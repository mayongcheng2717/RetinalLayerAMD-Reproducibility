#!/usr/bin/env python3
"""Independent numerical QA for consolidated v6 results.

Reads only verified returned results and read-only uploaded inputs. Writes a JSON
summary and a check-level CSV under work_v6/consolidated_results/qa.
"""
from pathlib import Path
import hashlib, io, json, tarfile, math
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path('/mnt/agents/output/审稿修订_本地补分析包')
WORK = ROOT / 'work_v6'
CON = WORK / 'consolidated_results'
OUT = CON / 'qa'
OUT.mkdir(parents=True, exist_ok=True)
RES = ROOT / 'returned_results'
UP = Path('/mnt/agents/upload')
LAYERS = ['CSI', 'GCC', 'GCL', 'INL', 'IPL', 'OPL', 'PS', 'RNFL', 'RPE']
EARLY = ['INL', 'OPL', 'RPE']
NB = 10000
checks = []

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()

def add(group, name, passed, detail='', max_diff=None):
    checks.append({'group': group, 'check': name, 'passed': bool(passed),
                   'detail': str(detail), 'max_abs_difference': max_diff})

def maxdiff(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return 0.0
    return float(np.max(np.abs(a[ok] - b[ok])))

# Manifest and file integrity
manifest = json.loads((CON / 'consolidated_analysis_manifest.json').read_text(encoding='utf-8'))
hash_fail = []
for name, meta in manifest['files'].items():
    path = CON / name
    if not path.exists() or sha(path) != meta['sha256']:
        hash_fail.append(name)
add('integrity', 'manifest_sha256', not hash_fail, hash_fail or '15 output files matched')
for name, meta in manifest['files'].items():
    if 'rows' in meta:
        x = pd.read_csv(CON / name)
        add('integrity', f'dimensions:{name}', x.shape == (meta['rows'], meta['columns']), x.shape)

# SMR: read six raw result files again, not the consolidated table
frames = []
for disease in ['iamdgc_lateamd', 'finngen_WET_AMD', 'finngen_DRY_AMD',
                'finngen_DM_MACULOPATHY', 'finngen_H7_RETIVASCOCCLUSION']:
    z = pd.read_csv(UP / f'{disease}.smr', sep='\t')
    z.insert(0, 'disease', disease)
    frames.append(z)
with tarfile.open(UP / '机制分析补料2.tar.gz', 'r:gz') as tf:
    raw = tf.extractfile('续跑/smr_winkler/results/winkler_earlyamd.smr').read()
z = pd.read_csv(io.BytesIO(raw), sep='\t')
z.insert(0, 'disease', 'winkler_earlyamd')
frames.append(z)
smr_raw = pd.concat(frames, ignore_index=True)
smr_out = pd.read_csv(CON / 'SMR_global_correction_2354.csv')
keys = ['disease', 'probeID']
smr_raw = smr_raw.sort_values(keys).reset_index(drop=True)
smr_out = smr_out.sort_values(keys).reset_index(drop=True)
add('SMR', 'raw_rows_and_unique_keys', len(smr_raw) == 2354 and not smr_raw.duplicated(keys).any(), len(smr_raw))
add('SMR', 'raw_output_key_match', smr_raw[keys].equals(smr_out[keys]))
for col in ['p_SMR', 'b_SMR', 'se_SMR', 'p_HEIDI']:
    d = maxdiff(smr_raw[col], smr_out[col])
    add('SMR', f'raw_output:{col}', d < 1e-12, max_diff=d)
# Independent BH step-up implementation
p = smr_raw.p_SMR.to_numpy(float)
m = len(p)
order = np.argsort(p, kind='mergesort')
ranked = p[order] * m / np.arange(1, m + 1)
q_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
q = np.empty(m)
q[order] = np.minimum(q_sorted, 1.0)
d = maxdiff(q, smr_out.q_BH_global2354)
add('SMR', 'independent_BH_2354', d < 1e-12, max_diff=d)
d = maxdiff(np.minimum(p * m, 1), smr_out.p_Bonferroni_global2354)
add('SMR', 'independent_global_Bonferroni', d < 1e-12, max_diff=d)
n_by_disease = smr_raw.groupby('disease').size().to_dict()
within = np.array([p[i] < 0.05 / n_by_disease[d] for i, d in enumerate(smr_raw.disease)])
heidi01 = smr_raw.p_HEIDI.to_numpy() > 0.01
heidi05 = smr_raw.p_HEIDI.to_numpy() > 0.05
global_bonf = p * m < 0.05
global_bh = q < 0.05
legacy = within & heidi01
rules = {
    'legacy_within_disease_Bonferroni_HEIDI_0_01': legacy,
    'legacy42_with_HEIDI_0_05': legacy & heidi05,
    'global_Bonferroni2354_HEIDI_0_01': global_bonf & heidi01,
    'global_Bonferroni2354_HEIDI_0_05': global_bonf & heidi05,
    'global_BH2354_HEIDI_0_01': global_bh & heidi01,
    'global_BH2354_HEIDI_0_05': global_bh & heidi05,
}
expected = {'legacy_within_disease_Bonferroni_HEIDI_0_01': (42, 27),
            'legacy42_with_HEIDI_0_05': (24, 18),
            'global_Bonferroni2354_HEIDI_0_01': (29, 19),
            'global_Bonferroni2354_HEIDI_0_05': (15, 11),
            'global_BH2354_HEIDI_0_01': (98, 63),
            'global_BH2354_HEIDI_0_05': (67, 45)}
for rule, mask in rules.items():
    got = (int(mask.sum()), int(smr_raw.loc[mask, 'probeID'].nunique()))
    add('SMR', f'counts:{rule}', got == expected[rule], got)
add('SMR', 'HEIDI_0_01_to_0_05_in_legacy42', int((legacy & (smr_raw.p_HEIDI <= 0.05)).sum()) == 18)
old = pd.read_csv(ROOT / 'inputs' / 'S1_SMR_42pairs.csv')
old_keys = set(zip(old.disease, old.Gene))
new_keys = set(zip(smr_raw.loc[legacy, 'disease'], smr_raw.loc[legacy, 'probeID']))
add('SMR', 'legacy42_exact_key_reproduction', old_keys == new_keys, f'{len(old_keys)} vs {len(new_keys)}')

# MVMR independent WLS via QR least squares on sqrt-weight transformed matrices
mm = pd.read_csv(RES / '02_mvmr_ucsc_verified/matrices/winkler_earlyamd.tsv.gz', sep='\t')
BX = mm[[f'bx_{x}' for x in EARLY]].to_numpy(float)
BXSE = mm[[f'bxse_{x}' for x in EARLY]].to_numpy(float)
by = mm.by.to_numpy(float)
byse = mm.byse.to_numpy(float)

def fit_wls(X, y, se):
    w = 1 / se**2
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    resid = y - X @ beta
    q = float(np.sum(w * resid**2))
    df = len(y) - X.shape[1]
    phi = max(q / df, 1.0)
    xtx_inv = np.linalg.pinv(X.T @ (w[:, None] * X))
    cov = xtx_inv * phi
    ses = np.sqrt(np.diag(cov))
    pv = 2 * stats.t.sf(np.abs(beta / ses), df)
    return beta, ses, pv, q, df, phi, resid, xtx_inv, w

def cond_f(X, XSE):
    ans = []
    for j in range(X.shape[1]):
        others = [k for k in range(X.shape[1]) if k != j]
        # Independent OLS residualization using statsmodels-like projection through pinv
        Xo = X[:, others]
        fitted = Xo @ (np.linalg.pinv(Xo.T @ Xo) @ Xo.T @ X[:, j])
        r = X[:, j] - fitted
        delta = np.linalg.pinv(Xo.T @ Xo) @ Xo.T @ X[:, j]
        var = XSE[:, j]**2 + np.sum((XSE[:, others] * delta)**2, axis=1)
        ans.append(float(np.mean(r**2 / var)))
    return np.array(ans)

beta, se, pv, q, df, phi, resid, xtx_inv, w = fit_wls(BX, by, byse)
f = cond_f(BX, BXSE)
model = pd.read_csv(CON / 'MVMR_earlyAMD_validated_model.csv').set_index('layer').loc[EARLY]
for label, got, target in [('beta', beta, model.beta), ('se', se, model.se), ('p_t', pv, model.p_t),
                           ('conditional_F', f, model.conditional_F_gencov0)]:
    d = maxdiff(got, target)
    add('MVMR', f'full_model:{label}', d < 1e-9, max_diff=d)
add('MVMR', 'full_model:Q', abs(q - model.Q.iloc[0]) < 1e-9, q, abs(q - model.Q.iloc[0]))
crit = stats.t.ppf(0.975, df)
add('MVMR', 'full_model:t_CI', max(maxdiff(beta-crit*se, model.ci95_low), maxdiff(beta+crit*se, model.ci95_high)) < 1e-9)
ref = pd.read_csv(RES / '02_mvmr_ucsc_verified/S7_mvmr_complete.csv')
ref = ref[ref.outcome.eq('winkler_earlyamd')].set_index('layer').loc[EARLY]
add('MVMR', 'three_level_reference_beta', maxdiff(beta, ref.mvmr_beta) < 1e-6, max_diff=maxdiff(beta, ref.mvmr_beta))
add('MVMR', 'three_level_reference_conditional_F', maxdiff(f, ref.conditional_F_gencov0) < 1e-6, max_diff=maxdiff(f, ref.conditional_F_gencov0))
# Full LOO grid
loo = pd.read_csv(CON / 'MVMR_earlyAMD_LOO.csv')
calc = []
for i in range(len(mm)):
    keep = np.arange(len(mm)) != i
    b, s, p0, q0, d0, ph0, r0, x0, w0 = fit_wls(BX[keep], by[keep], byse[keep])
    f0 = cond_f(BX[keep], BXSE[keep])
    for j, layer in enumerate(EARLY):
        calc.append((mm.SNP.iloc[i], layer, b[j], s[j], p0[j], q0, f0[j]))
calc = pd.DataFrame(calc, columns=['omitted_SNP', 'layer', 'b', 's', 'p', 'q', 'f'])
joined = loo.merge(calc, on=['omitted_SNP', 'layer'], how='outer', indicator=True)
add('MVMR', 'LOO_354_key_match', len(joined) == 354 and joined._merge.eq('both').all(), len(joined))
for a, b in [('beta', 'b'), ('se', 's'), ('p_t', 'p'), ('Q', 'q'), ('conditional_F_gencov0', 'f')]:
    d = maxdiff(joined[a], joined[b])
    add('MVMR', f'LOO:{a}', d < 1e-8, max_diff=d)
# Influence and collinearity
leverage = np.einsum('ij,jk,ik->i', BX, xtx_inv, BX) * w
qcon = w * resid**2
student = np.sqrt(w) * resid / np.sqrt(1-leverage)
cook = qcon * leverage / (3 * phi * (1-leverage)**2)
inf = pd.read_csv(CON / 'MVMR_earlyAMD_influence.csv').set_index('SNP').loc[mm.SNP]
for name, got, target in [('Q_contribution', qcon, inf.Q_contribution), ('leverage', leverage, inf.leverage),
                          ('studentized_residual', student, inf.studentized_residual), ('Cook_distance', cook, inf.Cook_distance)]:
    d = maxdiff(got, target)
    add('MVMR', f'influence:{name}', d < 1e-8, max_diff=d)
xtx = BX.T @ (w[:, None] * BX)
D = np.diag(1 / np.sqrt(np.diag(xtx)))
cor = D @ xtx @ D
vif = np.diag(np.linalg.pinv(cor))
eig = np.linalg.eigvalsh(cor)
coll = pd.read_csv(CON / 'MVMR_earlyAMD_collinearity.csv').set_index('layer').loc[EARLY]
add('MVMR', 'weighted_VIF', maxdiff(vif, coll.weighted_VIF) < 1e-9, max_diff=maxdiff(vif, coll.weighted_VIF))
add('MVMR', 'maximum_condition_index', abs(np.sqrt(eig.max()/eig.min()) - coll.maximum_condition_index.iloc[0]) < 1e-9)
rpe_key = loo[(loo.layer.eq('RPE')) & (loo.omitted_SNP.eq('rs10737680'))].iloc[0]
add('MVMR', 'RPE_rs10737680_dependence', abs(rpe_key.beta + 0.049356178409909156) < 1e-6 and abs(rpe_key.p_t - 0.5905612350331819) < 1e-6,
    f'beta={rpe_key.beta}; P={rpe_key.p_t}')
other = loo[(loo.layer.eq('RPE')) & (~loo.omitted_SNP.eq('rs10737680'))]
add('MVMR', 'RPE_other_117_LOO', len(other) == 117 and (other.p_t < 0.05).all(), int((other.p_t < 0.05).sum()))
idx = mm.index[mm.SNP.eq('rs10737680')][0]
add('MVMR', 'rs10737680_Q_fraction', abs(100*qcon[idx]/q - 36.730763632) < 1e-6, 100*qcon[idx]/q)

# Sample overlap scenarios
s3 = pd.read_csv(RES / '09_full54_scdrs_complete/results_for_return/MR/S3_full54_univariate_complete.csv')
early = s3[s3.outcome.eq('winkler_earlyamd')].copy()
sc = pd.read_csv(CON / 'earlyAMD_sample_overlap_bias_scenarios.csv')
expected_rows = len(early) * 2 * 5 * 4
add('overlap', 'scenario_dimensions', len(sc) == expected_rows == 360, len(sc))
max_bias_diff = max_shift_diff = 0.0
for _, r in sc.iterrows():
    source = early[early.layer.eq(r.layer)].iloc[0]
    f0 = source.F_mean if r.strength_parameter == 'mean_F' else source.F_min
    bias = r.absolute_observational_logOR_assumption * r.relative_overlap_fraction / f0
    sign = np.sign(source.IVW_beta)
    shifted = source.IVW_beta - sign*bias
    max_bias_diff = max(max_bias_diff, abs(bias-r.approximate_absolute_bias))
    max_shift_diff = max(max_shift_diff, abs(shifted-r.worst_toward_null_shifted_beta))
add('overlap', 'bias_formula_all_360', max_bias_diff < 1e-12, max_diff=max_bias_diff)
add('overlap', 'toward_null_shift_all_360', max_shift_diff < 1e-12, max_diff=max_shift_diff)
bounds = json.loads((CON / 'earlyAMD_sample_overlap_bounds.json').read_text(encoding='utf-8'))
accounting = (bounds['earlyAMD_cases'] + bounds['earlyAMD_controls'] == 105248 and
              bounds['earlyAMD_UKBB_cases'] + bounds['earlyAMD_UKBB_controls'] == 57802 and
              bounds['maximum_possible_shared_participants'] == 38485 and
              bounds['actual_shared_participants'] is None)
add('overlap', 'sample_accounting_and_unknown_actual_overlap', accounting)
extreme = sc[(sc.strength_parameter.eq('minimum_single_SNP_F')) &
             (sc.relative_overlap_fraction.eq(1.0)) &
             (sc.absolute_observational_logOR_assumption.eq(1.0)) & sc.layer.isin(['INL', 'OPL'])]
add('overlap', 'INL_OPL_extreme_scenario', len(extreme) == 2 and extreme.shifted_ci_excludes_zero.all(),
    extreme[['layer', 'worst_toward_null_shifted_beta', 'worst_toward_null_shifted_ci95_low', 'worst_toward_null_shifted_ci95_high']].to_dict('records'))

# Reverse MR: independent slope through origin and complete LOO/Wald checks
rev = pd.read_csv(CON / 'reverseMR_earlyAMD_unified.csv').set_index('layer')
revloo = pd.read_csv(CON / 'reverseMR_earlyAMD_LOO.csv')
wald = pd.read_csv(CON / 'reverseMR_earlyAMD_single_SNP_Wald.csv')
outtab = pd.read_csv(CON / 'reverseMR_earlyAMD_MR_PRESSO_outlier_tests.csv')

def slope_origin(x, y, yse):
    w = 1/yse**2
    b = float(np.sum(w*x*y)/np.sum(w*x*x))
    res = y-b*x
    q = float(np.sum(w*res**2))
    df = len(x)-1
    se = math.sqrt((q/df)/np.sum(w*x*x))
    return b, se, 2*stats.t.sf(abs(b/se), df), q, df

def primary_origin(x, y, yse):
    w = 1/yse**2
    b = float(np.sum(w*x*y)/np.sum(w*x*x))
    q = float(np.sum(w*(y-b*x)**2))
    df = len(x)-1
    se = math.sqrt(max(q/df, 1)/np.sum(w*x*x))
    return b, se, 2*stats.t.sf(abs(b/se), df), q, df

# Naive per-simulation MR-PRESSO recomputation, independent of vectorized generator
max_presso = {'global_RSS': 0.0, 'global_p': 0.0, 'raw_p': 0.0, 'corrected_beta': 0.0,
              'corrected_se': 0.0, 'corrected_p': 0.0, 'distortion': 0.0, 'distortion_p': 0.0}
with tarfile.open(UP / '机制分析补料2.tar.gz', 'r:gz') as tf:
    for lex_index, layer in enumerate(sorted(LAYERS)):
        member = f'本地续跑包/work/mr/harm/rev_winkler_earlyamd__{layer}.tsv'
        bbytes = tf.extractfile(member).read()
        dat = pd.read_csv(io.BytesIO(bbytes), sep='\t')
        bx0, by0 = dat.BETA_EXP.to_numpy(float), dat.BETA_OUT.to_numpy(float)
        bxse0, byse = dat.SE_EXP.to_numpy(float), dat.SE_OUT.to_numpy(float)
        orient = np.sign(bx0)
        bx, by, bxse = bx0*orient, by0*orient, bxse0
        n = len(bx)
        w = 1/byse**2
        r = rev.loc[layer]
        b, s, p0, q0, df0 = primary_origin(bx, by, byse)
        for label, got, target in [('beta', b, r.IVW_beta), ('se', s, r.IVW_se), ('p', p0, r.IVW_p_recalculated), ('Q', q0, r.IVW_Q)]:
            d = abs(got-target)
            add('reverseMR', f'{layer}:primary_{label}', d < 1e-8, max_diff=d)
        add('reverseMR', f'{layer}:input_sha256', hashlib.sha256(bbytes).hexdigest() == r.input_sha256)
        # Observed LOO slopes and residuals
        theta = np.empty(n)
        dif = np.empty(n)
        for j in range(n):
            keep = np.arange(n) != j
            theta[j] = np.sum(w[keep]*bx[keep]*by[keep])/np.sum(w[keep]*bx[keep]**2)
            dif[j] = by[j]-bx[j]*theta[j]
        rssobs = float(np.sum(w*dif**2))
        rng = np.random.default_rng(20261025 + 2000 + lex_index)
        gc = 0
        oc = np.zeros(n, dtype=int)
        completed = 0
        while completed < NB:
            bs = min(250, NB-completed)
            xr = rng.normal(bx[None, :], bxse[None, :], size=(bs, n))
            yr = rng.normal((bx*theta)[None, :], byse[None, :], size=(bs, n))
            for k in range(bs):
                rss = 0.0
                for j in range(n):
                    keep = np.arange(n) != j
                    slope = np.sum(w[keep]*xr[k, keep]*yr[k, keep])/np.sum(w[keep]*xr[k, keep]**2)
                    rr = yr[k, j]-xr[k, j]*slope
                    rss += w[j]*rr**2
                    expected_residual = yr[k, j]-xr[k, j]*theta[j]
                    oc[j] += int(expected_residual**2 > dif[j]**2)
                gc += int(rss > rssobs)
            completed += bs
        gp = gc/NB
        rawp = oc/NB
        adjp = np.minimum(rawp*n, 1)
        conduct = gp < 0.05
        outliers = np.flatnonzero(adjp <= 0.05) if conduct else np.array([], dtype=int)
        max_presso['global_RSS'] = max(max_presso['global_RSS'], abs(rssobs-r.MR_PRESSO_global_RSS_obs))
        max_presso['global_p'] = max(max_presso['global_p'], abs(gp-r.MR_PRESSO_global_p))
        expected_snps = ';'.join(dat.SNP.iloc[outliers].astype(str))
        got_snps = '' if pd.isna(r.MR_PRESSO_outlier_SNPs) else r.MR_PRESSO_outlier_SNPs
        add('reverseMR', f'{layer}:PRESSO_global_and_outliers', abs(gp-r.MR_PRESSO_global_p) < 1e-12 and expected_snps == got_snps,
            f'P={gp}; outliers={expected_snps}')
        raw_fit = slope_origin(bx, by, byse)
        corrected = (np.nan, np.nan, np.nan, np.nan, np.nan)
        dist = distp = np.nan
        if conduct and 0 < len(outliers) < n:
            keep = np.setdiff1d(np.arange(n), outliers)
            corrected = slope_origin(bx[keep], by[keep], byse[keep])
            biasobs = (raw_fit[0]-corrected[0])/abs(corrected[0])
            dist = 100*biasobs
            nonout = np.setdiff1d(np.arange(n), outliers)
            biase = np.empty(NB)
            for sim in range(NB):
                random_nonout = rng.choice(nonout, size=n-len(outliers), replace=True)
                inds = np.concatenate([outliers, random_nonout])[:n-len(outliers)]
                mod = slope_origin(bx[inds], by[inds], byse[inds])
                biase[sim] = (raw_fit[0]-mod[0])/abs(mod[0])
            distp = float(np.mean(np.abs(biase) > abs(biasobs)))
            for key, got, target in [('corrected_beta', corrected[0], r.MR_PRESSO_corrected_beta),
                                     ('corrected_se', corrected[1], r.MR_PRESSO_corrected_se),
                                     ('corrected_p', corrected[2], r.MR_PRESSO_corrected_p),
                                     ('distortion', dist, r.MR_PRESSO_distortion_coefficient_percent),
                                     ('distortion_p', distp, r.MR_PRESSO_distortion_p)]:
                max_presso[key] = max(max_presso[key], abs(got-target))
        if conduct:
            sub = outtab[outtab.layer.eq(layer)].set_index('SNP').loc[dat.SNP]
            max_presso['raw_p'] = max(max_presso['raw_p'], maxdiff(rawp, sub.Pvalue_raw))
            add('reverseMR', f'{layer}:PRESSO_outlier_table', len(sub) == n and maxdiff(adjp, sub.Pvalue_Bonferroni) < 1e-12)
        # LOO and Wald independent checks
        for j, snp in enumerate(dat.SNP):
            keep = np.arange(n) != j
            fit = primary_origin(bx[keep], by[keep], byse[keep])
            lr = revloo[(revloo.layer.eq(layer)) & (revloo.omitted_SNP.eq(snp))].iloc[0]
            if max(abs(fit[0]-lr.ivw_beta), abs(fit[1]-lr.ivw_se), abs(fit[2]-lr.ivw_p_t), abs(fit[3]-lr.Q)) >= 1e-8:
                add('reverseMR', f'{layer}:{snp}:LOO', False)
            wr = wald[(wald.layer.eq(layer)) & (wald.SNP.eq(snp))].iloc[0]
            ratio = by[j]/bx[j]
            ratio_se = abs(byse[j]/bx[j])
            if max(abs(ratio-wr.wald_ratio_beta), abs(ratio_se-wr.wald_ratio_se)) >= 1e-10:
                add('reverseMR', f'{layer}:{snp}:Wald', False)
add('reverseMR', 'all_72_LOO_and_Wald', not any((not c['passed']) and (':LOO' in c['check'] or ':Wald' in c['check']) for c in checks))
for key, value in max_presso.items():
    add('reverseMR', f'PRESSO_independent_recalculation:{key}', value < 1e-8, max_diff=value)
ps = rev.loc['PS']
add('reverseMR', 'PS_formal_outlier_rule', ps.MR_PRESSO_global_p == 0.0616 and not ps.MR_PRESSO_outlier_test_conducted and ps.MR_PRESSO_n_outliers_Bonferroni_0_05 == 0)

# Cross-result RPE univariable omission check from full54 LOO output
loo54_path = RES / '09_full54_scdrs_complete/results_for_return/MR/MR_leave_one_out_full54.csv'
loo54 = pd.read_csv(loo54_path)
rpe_full = early[early.layer.eq('RPE')].iloc[0]
rpe_loo = loo54[(loo54.layer.eq('RPE')) & (loo54.outcome.eq('winkler_earlyamd')) & (loo54.omitted_SNP.eq('rs10737680'))]
add('cross_result', 'RPE_uvmr_single_SNP_dependence', len(rpe_loo) == 1 and abs(rpe_loo.iloc[0].ivw_p - 0.387772) < 1e-5 and abs(rpe_loo.iloc[0].ivw_beta + 0.173169) < 1e-5,
    rpe_loo.to_dict('records'))

# Robust MVMR: independent QuantileRegressor and statsmodels WLS checks
# Restore the 118-instrument arrays after the reverse-MR loop used local names.
by = mm.by.to_numpy(float)
byse = mm.byse.to_numpy(float)
from sklearn.linear_model import QuantileRegressor
import statsmodels.api as sm
med = pd.read_csv(CON / 'MVMR_earlyAMD_median_sensitivity.csv')
boot = pd.read_csv(CON / 'MVMR_earlyAMD_median_bootstrap.csv.gz')
add('robustMVMR', 'bootstrap_60000_rows', len(boot) == 60000 and not boot.beta_bootstrap.isna().any(), len(boot))
for model_name, indices in [('all_118_instruments', np.arange(len(mm))),
                            ('without_rs10737680', np.flatnonzero(~mm.SNP.eq('rs10737680').to_numpy()))]:
    qr = QuantileRegressor(quantile=0.5, alpha=0, fit_intercept=False, solver='highs')
    qr.fit(BX[indices], by[indices], sample_weight=1/byse[indices]**2)
    for j, layer in enumerate(EARLY):
        row = med[(med.model.eq(model_name)) & (med.layer.eq(layer))].iloc[0]
        d = abs(qr.coef_[j]-row.beta)
        add('robustMVMR', f'{model_name}:{layer}:median_point', d < 1e-7, max_diff=d)
        vals = boot[(boot.model.eq(model_name)) & (boot.layer.eq(layer))].beta_bootstrap.to_numpy()
        sd = vals.std(ddof=1)
        add('robustMVMR', f'{model_name}:{layer}:bootstrap_SE', len(vals) == 10000 and abs(sd-row.se_parametric_bootstrap) < 1e-10,
            max_diff=abs(sd-row.se_parametric_bootstrap))
        add('robustMVMR', f'{model_name}:{layer}:bootstrap_percentiles',
            max(abs(np.quantile(vals, .025)-row.bootstrap_percentile_ci95_low),
                abs(np.quantile(vals, .975)-row.bootstrap_percentile_ci95_high)) < 1e-10)
# Independent MVMR-Egger WLS fit
orient = np.sign(BX[:, 0])
design = np.c_[np.ones(len(by)), BX*orient[:, None]]
egger_wls = sm.WLS(by*orient, design, weights=1/byse**2).fit()
eg = pd.read_csv(CON / 'MVMR_earlyAMD_Egger_sensitivity.csv').set_index('term').loc[['Intercept']+EARLY]
add('robustMVMR', 'Egger_beta', maxdiff(egger_wls.params, eg.beta) < 1e-8, max_diff=maxdiff(egger_wls.params, eg.beta))
add('robustMVMR', 'Egger_SE', maxdiff(egger_wls.bse, eg.se) < 1e-8, max_diff=maxdiff(egger_wls.bse, eg.se))
q_eg = float(np.sum((1/byse**2)*egger_wls.resid**2))
add('robustMVMR', 'Egger_Q', abs(q_eg-eg.Q.iloc[0]) < 1e-8, max_diff=abs(q_eg-eg.Q.iloc[0]))
add('robustMVMR', 'RPE_robust_direction_discordance',
    abs(med[(med.model.eq('all_118_instruments')) & med.layer.eq('RPE')].beta.iloc[0]-0.1238959) < 1e-6 and
    med[(med.model.eq('all_118_instruments')) & med.layer.eq('RPE')].p_normal.iloc[0] > 0.05)

frame = pd.DataFrame(checks)
# Collapse absent explicit per-SNP failure records into the group checks; no synthetic passes added.
failed = frame[~frame.passed]
summary = {'status': 'PASS' if failed.empty else 'FAIL', 'checks': len(frame),
           'passed': int(frame.passed.sum()), 'failed': int((~frame.passed).sum()),
           'failed_checks': failed.to_dict('records'),
           'scope': {'SMR_tests': 2354, 'MVMR_LOO_rows': 354, 'overlap_scenarios': 360,
                     'reverse_MR_models': 9, 'reverse_LOO_rows': 72, 'reverse_Wald_rows': 72,
                     'MVMR_median_bootstrap_rows': 60000, 'MVMR_Egger_terms': 4},
           'independence': 'Recomputed from original uploaded and returned inputs using separate QR/naive-loop implementations; generator assertions were not used as evidence.'}
frame.to_csv(OUT / '整合分析_独立验收_逐项检查.csv', index=False)
(OUT / '整合分析_独立验收报告.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(summary, ensure_ascii=False, indent=2))
if not failed.empty:
    raise SystemExit(1)
