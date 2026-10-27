#!/usr/bin/env python3
"""Consolidated v6 analyses from verified returned data.

This script does not download data or alter any uploaded file. It writes only to
work_v6/consolidated_results. All numerical outputs are tied to fixed inputs and
independent validation checks.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path('/mnt/agents/output/审稿修订_本地补分析包')
UPLOAD = Path('/mnt/agents/upload')
OUT = ROOT / 'work_v6' / 'consolidated_results'
OUT.mkdir(parents=True, exist_ok=True)
RES = ROOT / 'returned_results'
LAYERS = ['CSI', 'GCC', 'GCL', 'INL', 'IPL', 'OPL', 'PS', 'RNFL', 'RPE']
EARLY_LAYERS = ['INL', 'OPL', 'RPE']
BASE_SEED = 20261025
NB = 10000
BATCH = 250
ALPHA = 0.05
TOL = 1e-6
manifest = {'analyses': {}, 'validation': {}, 'files': {}}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save_csv(frame: pd.DataFrame, name: str) -> None:
    path = OUT / name
    frame.to_csv(path, index=False)
    manifest['files'][name] = {'rows': len(frame), 'columns': len(frame.columns),
                               'sha256': sha256_file(path)}


def save_json(obj: dict, name: str) -> None:
    path = OUT / name
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')
    manifest['files'][name] = {'sha256': sha256_file(path)}


# 1. Complete SMR second-stage global correction
smr_frames = []
for disease in ['iamdgc_lateamd', 'finngen_WET_AMD', 'finngen_DRY_AMD',
                'finngen_DM_MACULOPATHY', 'finngen_H7_RETIVASCOCCLUSION']:
    frame = pd.read_csv(UPLOAD / f'{disease}.smr', sep='\t')
    frame.insert(0, 'disease', disease)
    smr_frames.append(frame)
with tarfile.open(UPLOAD / '机制分析补料2.tar.gz', 'r:gz') as archive:
    member = '续跑/smr_winkler/results/winkler_earlyamd.smr'
    early_bytes = archive.extractfile(member).read()
early = pd.read_csv(io.BytesIO(early_bytes), sep='\t')
early.insert(0, 'disease', 'winkler_earlyamd')
smr_frames.append(early)
smr = pd.concat(smr_frames, ignore_index=True)
assert len(smr) == 2354
smr['n_tests_within_disease'] = smr.groupby('disease')['probeID'].transform('size')
smr['pass_within_bonferroni'] = smr.p_SMR < 0.05 / smr.n_tests_within_disease
smr['pass_HEIDI_0_01'] = smr.p_HEIDI > 0.01
smr['pass_HEIDI_0_05'] = smr.p_HEIDI > 0.05
smr['legacy42_selected'] = smr.pass_within_bonferroni & smr.pass_HEIDI_0_01
smr['p_Bonferroni_global2354'] = np.minimum(1.0, smr.p_SMR * len(smr))
smr['q_BH_global2354'] = multipletests(smr.p_SMR, method='fdr_bh')[1]
smr['pass_global_Bonferroni'] = smr.p_Bonferroni_global2354 < 0.05
smr['pass_global_BH'] = smr.q_BH_global2354 < 0.05
legacy = pd.read_csv(ROOT / 'inputs' / 'S1_SMR_42pairs.csv')
legacy_keys = set(zip(legacy.disease, legacy.Gene))
calculated_keys = set(zip(smr.loc[smr.legacy42_selected, 'disease'],
                          smr.loc[smr.legacy42_selected, 'probeID']))
assert calculated_keys == legacy_keys
legacy_check = legacy.merge(
    smr.loc[smr.legacy42_selected, ['disease', 'probeID', 'p_SMR', 'p_HEIDI']],
    left_on=['disease', 'Gene'], right_on=['disease', 'probeID'], how='outer', indicator=True)
assert legacy_check['_merge'].eq('both').all()
assert np.max(np.abs(legacy_check.p_SMR_x - legacy_check.p_SMR_y)) == 0
assert np.max(np.abs(legacy_check.p_HEIDI_x - legacy_check.p_HEIDI_y)) == 0
symbol_map = legacy.drop_duplicates('Gene').set_index('Gene')['symbol'].to_dict()
smr['symbol_legacy42_annotation'] = smr.probeID.map(symbol_map)
smr['second_stage_scope'] = (
    'Global correction covers 2354 disease-probe tests entering SMR Select; '
    'it does not correct the preceding region-selection stage.')
save_csv(smr, 'SMR_global_correction_2354.csv')
rules = {
    'legacy_within_disease_Bonferroni_HEIDI_0_01': smr.legacy42_selected,
    'legacy42_with_HEIDI_0_05': smr.legacy42_selected & smr.pass_HEIDI_0_05,
    'global_Bonferroni2354_HEIDI_0_01': smr.pass_global_Bonferroni & smr.pass_HEIDI_0_01,
    'global_Bonferroni2354_HEIDI_0_05': smr.pass_global_Bonferroni & smr.pass_HEIDI_0_05,
    'global_BH2354_HEIDI_0_01': smr.pass_global_BH & smr.pass_HEIDI_0_01,
    'global_BH2354_HEIDI_0_05': smr.pass_global_BH & smr.pass_HEIDI_0_05,
}
count_rows = []
for rule, mask in rules.items():
    selected = smr.loc[mask]
    old = selected.legacy42_selected
    count_rows.append({
        'rule': rule, 'pairs': len(selected), 'unique_probes': selected.probeID.nunique(),
        'pairs_from_legacy42': int(old.sum()), 'new_pairs_vs_legacy42': int((~old).sum()),
        'unique_new_probes_vs_legacy42': selected.loc[~old, 'probeID'].nunique(),
        'interpretation': ('Legacy candidate catalogue' if rule.startswith('legacy') else
                           'Second-stage global sensitivity; BH additions require separate first-stage and locus review')})
counts = pd.DataFrame(count_rows)
save_csv(counts, 'SMR_global_correction_counts.csv')
expected = {
    'legacy_within_disease_Bonferroni_HEIDI_0_01': (42, 27),
    'legacy42_with_HEIDI_0_05': (24, None),
    'global_Bonferroni2354_HEIDI_0_01': (29, 19),
    'global_Bonferroni2354_HEIDI_0_05': (15, 11),
    'global_BH2354_HEIDI_0_01': (98, 63),
    'global_BH2354_HEIDI_0_05': (67, 45),
}
for row in count_rows:
    ep, eg = expected[row['rule']]
    assert row['pairs'] == ep
    if eg is not None:
        assert row['unique_probes'] == eg
legacy_global = smr.loc[smr.legacy42_selected].copy()
legacy_global['global_Bonferroni_retained_HEIDI_0_01'] = (
    legacy_global.pass_global_Bonferroni & legacy_global.pass_HEIDI_0_01)
legacy_global['global_Bonferroni_retained_HEIDI_0_05'] = (
    legacy_global.pass_global_Bonferroni & legacy_global.pass_HEIDI_0_05)
legacy_global['HEIDI_0_01_to_0_05'] = (
    (legacy_global.p_HEIDI > 0.01) & (legacy_global.p_HEIDI <= 0.05))
assert legacy_global.HEIDI_0_01_to_0_05.sum() == 18
save_csv(legacy_global, 'SMR_legacy42_global_sensitivity.csv')
manifest['analyses']['SMR'] = {
    'tests': 2354, 'legacy42_exactly_reproduced': True,
    'global_Bonferroni_HEIDI_0_01': {'pairs': 29, 'unique_probes': 19},
    'global_Bonferroni_HEIDI_0_05': {'pairs': 15, 'unique_probes': 11},
    'global_BH_HEIDI_0_01': {'pairs': 98, 'unique_probes': 63, 'new_pairs': 56},
    'HEIDI_0_01_to_0_05_within_legacy42': 18,
    'two_stage_limitation_retained': True}


# 2. Early-AMD MVMR model validation, collinearity, LOO and influence
matrix_path = RES / '02_mvmr_ucsc_verified' / 'matrices' / 'winkler_earlyamd.tsv.gz'
mm = pd.read_csv(matrix_path, sep='\t')
reference = pd.read_csv(RES / '02_mvmr_ucsc_verified' / 'S7_mvmr_complete.csv')
reference = reference[reference.outcome.eq('winkler_earlyamd')].set_index('layer').loc[EARLY_LAYERS]
BX = mm[[f'bx_{x}' for x in EARLY_LAYERS]].to_numpy(float)
BXSE = mm[[f'bxse_{x}' for x in EARLY_LAYERS]].to_numpy(float)
by = mm.by.to_numpy(float)
byse = mm.byse.to_numpy(float)


def fit_mvmr(BX, by, byse):
    weights = 1 / byse ** 2
    xtx = BX.T @ (weights[:, None] * BX)
    covariance0 = np.linalg.inv(xtx)
    beta = covariance0 @ (BX.T @ (weights * by))
    residual = by - BX @ beta
    q = float(np.sum(weights * residual ** 2))
    df = len(by) - BX.shape[1]
    phi = max(q / df, 1.0)
    se = np.sqrt(np.diag(covariance0 * phi))
    p = 2 * stats.t.sf(np.abs(beta / se), df)
    return {'beta': beta, 'se': se, 'p': p, 'q': q, 'df': df, 'phi': phi,
            'q_p': float(stats.chi2.sf(q, df)), 'residual': residual,
            'covariance0': covariance0, 'weights': weights}


def conditional_f_gencov0(BX, BXSE):
    output = []
    for j in range(BX.shape[1]):
        other = [k for k in range(BX.shape[1]) if k != j]
        delta = np.linalg.lstsq(BX[:, other], BX[:, j], rcond=None)[0]
        residual = BX[:, j] - BX[:, other] @ delta
        sigma2 = BXSE[:, j] ** 2 + np.sum(BXSE[:, other] ** 2 * delta ** 2, axis=1)
        output.append(float(np.mean(residual ** 2 / sigma2)))
    return np.asarray(output)


full = fit_mvmr(BX, by, byse)
full_f = conditional_f_gencov0(BX, BXSE)
for j, layer in enumerate(EARLY_LAYERS):
    ref = reference.loc[layer]
    for value, target in [(full['beta'][j], ref.mvmr_beta), (full['se'][j], ref.mvmr_se),
                          (full['p'][j], ref.mvmr_p_t), (full_f[j], ref.conditional_F_gencov0),
                          (full['q'], ref.q), (full['q_p'], ref.q_p)]:
        assert abs(float(value) - float(target)) <= TOL
model_rows = []
crit = stats.t.ppf(0.975, full['df'])
for j, layer in enumerate(EARLY_LAYERS):
    model_rows.append({'outcome': 'winkler_earlyamd', 'layer': layer, 'nsnp': len(mm),
                       'exposures': '|'.join(EARLY_LAYERS), 'beta': full['beta'][j],
                       'se': full['se'][j], 'df': full['df'], 'p_t': full['p'][j],
                       'ci95_low': full['beta'][j] - crit * full['se'][j],
                       'ci95_high': full['beta'][j] + crit * full['se'][j],
                       'conditional_F_gencov0': full_f[j], 'Q': full['q'],
                       'Q_df': full['df'], 'Q_p': full['q_p'],
                       'covariance_assumption': 'gencov=0; same-sample OCT exposure covariance is not modelled'})
save_csv(pd.DataFrame(model_rows), 'MVMR_earlyAMD_validated_model.csv')
weights = full['weights']
xtx = BX.T @ (weights[:, None] * BX)
scale = np.diag(1 / np.sqrt(np.diag(xtx)))
cor = scale @ xtx @ scale
eigenvalues = np.linalg.eigvalsh(cor)
vif = np.diag(np.linalg.inv(cor))
collinearity = pd.DataFrame({'layer': EARLY_LAYERS, 'weighted_VIF': vif})
collinearity['weighted_condition_number'] = eigenvalues.max() / eigenvalues.min()
collinearity['maximum_condition_index'] = math.sqrt(eigenvalues.max() / eigenvalues.min())
collinearity['interpretation'] = 'Design collinearity diagnostic, not a test of MR validity or pleiotropy'
save_csv(collinearity, 'MVMR_earlyAMD_collinearity.csv')
save_csv(pd.DataFrame({'eigenvalue': eigenvalues,
                       'condition_index': np.sqrt(eigenvalues.max() / eigenvalues)}),
         'MVMR_earlyAMD_weighted_design_eigenvalues.csv')
# SNP-effect correlations are descriptive and distinct from genome-wide LDSC rg.
ldsc = pd.read_csv(ROOT / 'inputs' / 'ldsc_rg.csv')
correlation_rows = []
for a in range(3):
    for b in range(a + 1, 3):
        la, lb = EARLY_LAYERS[a], EARLY_LAYERS[b]
        pearson = stats.pearsonr(BX[:, a], BX[:, b])
        spearman = stats.spearmanr(BX[:, a], BX[:, b])
        match = ldsc[((ldsc.trait1.eq(f'zekavat_{la}')) & (ldsc.trait2.eq(f'zekavat_{lb}'))) |
                     ((ldsc.trait2.eq(f'zekavat_{la}')) & (ldsc.trait1.eq(f'zekavat_{lb}')))]
        assert len(match) == 1
        rg = match.iloc[0]
        correlation_rows.append({'layer1': la, 'layer2': lb, 'n_instruments': len(mm),
                                 'instrument_effect_pearson_r': pearson.statistic,
                                 'instrument_effect_pearson_p': pearson.pvalue,
                                 'instrument_effect_spearman_rho': spearman.statistic,
                                 'instrument_effect_spearman_p': spearman.pvalue,
                                 'LDSC_rg': rg.rg, 'LDSC_rg_se': rg.rg_se, 'LDSC_p': rg.p,
                                 'note': 'Instrument-effect correlations and genome-wide LDSC rg are different quantities'})
save_csv(pd.DataFrame(correlation_rows), 'MVMR_earlyAMD_exposure_correlations.csv')
loo_rows = []
for i in range(len(mm)):
    keep = np.arange(len(mm)) != i
    fit = fit_mvmr(BX[keep], by[keep], byse[keep])
    conditional_f = conditional_f_gencov0(BX[keep], BXSE[keep])
    critical = stats.t.ppf(0.975, fit['df'])
    for j, layer in enumerate(EARLY_LAYERS):
        loo_rows.append({'outcome': 'winkler_earlyamd', 'omitted_SNP': mm.SNP.iloc[i],
                         'layer': layer, 'nsnp_remaining': len(mm) - 1,
                         'beta': fit['beta'][j], 'se': fit['se'][j], 'p_t': fit['p'][j],
                         'ci95_low': fit['beta'][j] - critical * fit['se'][j],
                         'ci95_high': fit['beta'][j] + critical * fit['se'][j],
                         'beta_change_from_full': fit['beta'][j] - full['beta'][j],
                         'abs_beta_change_in_full_SE': abs(fit['beta'][j] - full['beta'][j]) / full['se'][j],
                         'direction_matches_full': bool(np.sign(fit['beta'][j]) == np.sign(full['beta'][j])),
                         'conditional_F_gencov0': conditional_f[j], 'Q': fit['q'],
                         'Q_df': fit['df'], 'Q_p': fit['q_p']})
mvmr_loo = pd.DataFrame(loo_rows)
save_csv(mvmr_loo, 'MVMR_earlyAMD_LOO.csv')
leverage = np.einsum('ij,jk,ik->i', BX, full['covariance0'], BX) * weights
q_contribution = weights * full['residual'] ** 2
studentized = np.sqrt(weights) * full['residual'] / np.sqrt(1 - leverage)
cook = q_contribution * leverage / (3 * full['phi'] * (1 - leverage) ** 2)
influence = pd.DataFrame({'SNP': mm.SNP, 'CHR': mm.CHR, 'BP': mm.BP,
                          'Q_contribution': q_contribution,
                          'Q_contribution_percent': 100 * q_contribution / full['q'],
                          'leverage': leverage, 'studentized_residual': studentized,
                          'Cook_distance': cook})
max_shift = mvmr_loo.groupby('omitted_SNP').abs_beta_change_in_full_SE.max()
influence['max_abs_beta_change_in_full_SE'] = influence.SNP.map(max_shift)
influence = influence.sort_values('Cook_distance', ascending=False)
save_csv(influence, 'MVMR_earlyAMD_influence.csv')
key = mvmr_loo[(mvmr_loo.omitted_SNP.eq('rs10737680')) & (mvmr_loo.layer.eq('RPE'))].iloc[0]
assert abs(key.beta - (-0.049356178409909156)) < TOL
assert abs(key.p_t - 0.5905612350331819) < TOL
assert abs(float(influence.loc[influence.SNP.eq('rs10737680'), 'Q_contribution_percent'].iloc[0]) - 36.730763632) < 1e-6
manifest['analyses']['MVMR_earlyAMD'] = {
    'reference_exactly_reproduced_within_1e_6': True, 'instruments': 118,
    'conditional_F': dict(zip(EARLY_LAYERS, full_f.tolist())),
    'Q': full['q'], 'Q_df': full['df'], 'Q_p': full['q_p'],
    'weighted_VIF': dict(zip(EARLY_LAYERS, vif.tolist())),
    'rs10737680_Q_contribution_percent': 36.730763632,
    'RPE_without_rs10737680': {'beta': float(key.beta), 'se': float(key.se), 'p_t': float(key.p_t)},
    'RPE_other_LOO_significant_count': int(((mvmr_loo.layer.eq('RPE')) &
                                            (~mvmr_loo.omitted_SNP.eq('rs10737680')) &
                                            (mvmr_loo.p_t < 0.05)).sum()),
    'interpretation': 'RPE association is single-SNP dependent and residual heterogeneity remains after omission'}


# 3. Early-AMD sample-overlap bounds and approximate bias scenarios
s3 = pd.read_csv(RES / '09_full54_scdrs_complete' / 'results_for_return' / 'MR' /
                 'S3_full54_univariate_complete.csv')
early_mr = s3[s3.outcome.eq('winkler_earlyamd')].copy()
assert len(early_mr) == 9
n_oct = 38485
n_early_cases, n_early_controls = 14034, 91214
n_early = n_early_cases + n_early_controls
n_early_ukb_cases, n_early_ukb_controls = 3105, 54697
n_early_ukb = n_early_ukb_cases + n_early_ukb_controls
bounds = {
    'OCT_reported_N': n_oct, 'OCT_N_rule': 'median across variant-specific N values',
    'earlyAMD_cases': n_early_cases, 'earlyAMD_controls': n_early_controls,
    'earlyAMD_total_N': n_early, 'earlyAMD_UKBB_cases': n_early_ukb_cases,
    'earlyAMD_UKBB_controls': n_early_ukb_controls, 'earlyAMD_UKBB_total_N': n_early_ukb,
    'maximum_possible_shared_participants': min(n_oct, n_early_ukb),
    'maximum_exposure_side_overlap_fraction': min(n_oct, n_early_ukb) / n_oct,
    'maximum_outcome_side_overlap_fraction': min(n_oct, n_early_ukb) / n_early,
    'UKBB_fraction_of_earlyAMD_total': n_early_ukb / n_early,
    'actual_shared_participants': None,
    'sources': {
        'OCT_N': 'returned_results/06_sample_sizes/OCT_GWAS_sample_sizes_by_layer.csv',
        'earlyAMD_counts': 'Winkler et al. 2020, PMC7449002, Supplementary Table 1',
        'earlyAMD_UKBB_counts': 'Winkler et al. 2020, PMC7449002, Supplementary Table 1'},
    'caution': 'Participant maxima are accounting bounds, not estimates of actual overlap or an effective covariance coefficient.'}
save_json(bounds, 'earlyAMD_sample_overlap_bounds.json')
scenario_rows = []
for _, row in early_mr.iterrows():
    for strength_name, strength in [('mean_F', row.F_mean), ('minimum_single_SNP_F', row.F_min)]:
        for overlap in [0, 0.25, 0.50, 0.75, 1.0]:
            for observational in [0.1, 0.2, 0.5, 1.0]:
                bias = observational * overlap / strength
                sign = np.sign(row.IVW_beta)
                shifted_beta = row.IVW_beta - sign * bias
                shifted_low = row.IVW_ci_low - sign * bias
                shifted_high = row.IVW_ci_high - sign * bias
                scenario_rows.append({
                    'outcome': row.outcome, 'layer': row.layer, 'nsnp': row.nsnp,
                    'IVW_beta': row.IVW_beta, 'IVW_se': row.IVW_se,
                    'IVW_ci95_low': row.IVW_ci_low, 'IVW_ci95_high': row.IVW_ci_high,
                    'strength_parameter': strength_name, 'F_used': strength,
                    'relative_overlap_fraction': overlap,
                    'absolute_observational_logOR_assumption': observational,
                    'approximate_absolute_bias': bias,
                    'bias_percent_of_abs_IVW_beta': 100 * bias / abs(row.IVW_beta),
                    'worst_toward_null_shifted_beta': shifted_beta,
                    'worst_toward_null_shifted_ci95_low': shifted_low,
                    'worst_toward_null_shifted_ci95_high': shifted_high,
                    'shifted_ci_excludes_zero': bool((shifted_low > 0) or (shifted_high < 0)),
                    'formula': 'absolute bias = absolute observational logOR x overlap fraction / F',
                    'interpretation': 'Approximate Burgess-Davies-Thompson sensitivity, not an exact correction for binary outcome or unequal sample sizes'})
scenarios = pd.DataFrame(scenario_rows)
save_csv(scenarios, 'earlyAMD_sample_overlap_bias_scenarios.csv')
extreme = scenarios[(scenarios.strength_parameter.eq('minimum_single_SNP_F')) &
                    (scenarios.relative_overlap_fraction.eq(1.0)) &
                    (scenarios.absolute_observational_logOR_assumption.eq(1.0)) &
                    (scenarios.layer.isin(['INL', 'OPL']))]
assert extreme.shifted_ci_excludes_zero.all()
manifest['analyses']['sample_overlap'] = {
    **bounds, 'scenario_rows': len(scenarios),
    'extreme_scenario': 'minimum single-SNP F, complete relative overlap, absolute observational logOR=1',
    'INL_OPL_extreme_shifted_CI_excludes_zero': True,
    'extreme_results': extreme[['layer', 'approximate_absolute_bias',
                                'worst_toward_null_shifted_beta',
                                'worst_toward_null_shifted_ci95_low',
                                'worst_toward_null_shifted_ci95_high']].to_dict('records'),
    'non_overlapping_replication_available': False,
    'method_reference': 'Burgess, Davies and Thompson 2016, PMC5082560, DOI 10.1002/gepi.21998'}


# 4. Reverse early-AMD MR-PRESSO and influence using one harmonized input per model

def loo_statistics(bx, by, weights):
    sxx = np.sum(weights * bx * bx, axis=-1)
    sxy = np.sum(weights * bx * by, axis=-1)
    theta = (sxy[..., None] - weights * bx * by) / (sxx[..., None] - weights * bx * bx)
    difference = by - bx * theta
    rss = np.sum(weights * difference * difference, axis=-1)
    return theta, difference, rss


def weighted_lm_official(bx, by, byse):
    bx, by, byse = [np.asarray(x, float) for x in (bx, by, byse)]
    weights = 1 / byse ** 2
    theta = np.sum(weights * bx * by) / np.sum(weights * bx * bx)
    residual = by - theta * bx
    q = np.sum(weights * residual ** 2)
    df = len(bx) - 1
    se = math.sqrt((q / df) / np.sum(weights * bx * bx))
    return {'beta': theta, 'se': se, 'p': 2 * stats.t.sf(abs(theta / se), df),
            'q': q, 'df': df}


def ivw_primary(bx, by, byse):
    bx, by, byse = [np.asarray(x, float) for x in (bx, by, byse)]
    weights = 1 / byse ** 2
    theta = np.sum(weights * bx * by) / np.sum(weights * bx * bx)
    q = np.sum(weights * (by - theta * bx) ** 2)
    df = len(bx) - 1
    se = math.sqrt(max(q / df, 1) / np.sum(weights * bx * bx))
    return {'beta': theta, 'se': se, 'p': 2 * stats.t.sf(abs(theta / se), df),
            'q': q, 'df': df, 'q_p': float(stats.chi2.sf(q, df))}


def simulate_presso(bx, by, bxse, byse, seed):
    bx, by, bxse, byse = [np.asarray(x, float) for x in (bx, by, bxse, byse)]
    orient = np.sign(bx)
    by, bx = by * orient, bx * orient
    n = len(bx)
    if n <= 3 or n >= NB:
        raise ValueError('Invalid MR-PRESSO instrument count')
    weights = 1 / byse ** 2
    theta_loo, dif_obs, rss_obs = loo_statistics(bx, by, weights)
    rss_obs = float(rss_obs)
    rng = np.random.default_rng(seed)
    global_count = 0
    outlier_count = np.zeros(n, dtype=np.int64)
    completed = 0
    while completed < NB:
        batch = min(BATCH, NB - completed)
        bx_random = rng.normal(bx[None, :], bxse[None, :], size=(batch, n))
        by_random = rng.normal((bx * theta_loo)[None, :], byse[None, :], size=(batch, n))
        _, _, rss_random = loo_statistics(bx_random, by_random, weights)
        global_count += int(np.sum(rss_random > rss_obs))
        expected = by_random - bx_random * theta_loo[None, :]
        outlier_count += np.sum(expected ** 2 > dif_obs[None, :] ** 2, axis=0)
        completed += batch
    global_p = global_count / NB
    raw_p = outlier_count / NB
    adjusted_p = np.minimum(raw_p * n, 1)
    conduct = global_p < ALPHA
    outlier_idx = np.flatnonzero(adjusted_p <= ALPHA) if conduct else np.array([], dtype=int)
    raw = weighted_lm_official(bx, by, byse)
    corrected = {'beta': np.nan, 'se': np.nan, 'p': np.nan, 'q': np.nan, 'df': np.nan}
    distortion_coefficient = np.nan
    distortion_p = np.nan
    if conduct and 0 < len(outlier_idx) < n:
        keep = np.setdiff1d(np.arange(n), outlier_idx)
        corrected = weighted_lm_official(bx[keep], by[keep], byse[keep])
        bias_observed = (raw['beta'] - corrected['beta']) / abs(corrected['beta'])
        distortion_coefficient = 100 * bias_observed
        nonoutlier = np.setdiff1d(np.arange(n), outlier_idx)
        bias_expected = np.empty(NB)
        for simulation in range(NB):
            random_nonoutlier = rng.choice(nonoutlier, size=n - len(outlier_idx), replace=True)
            indices = np.concatenate([outlier_idx, random_nonoutlier])
            indices = indices[:n - len(outlier_idx)]
            model = weighted_lm_official(bx[indices], by[indices], byse[indices])
            bias_expected[simulation] = (raw['beta'] - model['beta']) / abs(model['beta'])
        distortion_p = float(np.mean(np.abs(bias_expected) > abs(bias_observed)))
    return {'global_RSS_obs': rss_obs, 'global_p': global_p,
            'global_p_report': f'<{1 / NB:g}' if global_p == 0 else str(global_p),
            'outlier_test_conducted': conduct, 'outlier_raw_p': raw_p,
            'outlier_adjusted_p': adjusted_p, 'outlier_idx': outlier_idx,
            'dif_obs_squared': dif_obs ** 2, 'raw': raw, 'corrected': corrected,
            'distortion_coefficient_percent': distortion_coefficient,
            'distortion_p': distortion_p}


reverse_reference = pd.read_csv(RES / '01_mr' / 'S4_reverseMR_earlyAMD_recalculated.csv')
reverse_rows, outlier_rows, reverse_loo_rows, wald_rows = [], [], [], []
with tarfile.open(UPLOAD / '机制分析补料2.tar.gz', 'r:gz') as archive:
    for index, layer in enumerate(sorted(LAYERS)):
        member = f'本地续跑包/work/mr/harm/rev_winkler_earlyamd__{layer}.tsv'
        raw_bytes = archive.extractfile(member).read()
        data = pd.read_csv(io.BytesIO(raw_bytes), sep='\t')
        ref = reverse_reference[reverse_reference.layer.eq(layer)]
        assert len(ref) == 1 and len(data) == int(ref.nsnp.iloc[0])
        ref = ref.iloc[0]
        bx, bxse = data.BETA_EXP.to_numpy(float), data.SE_EXP.to_numpy(float)
        by, byse = data.BETA_OUT.to_numpy(float), data.SE_OUT.to_numpy(float)
        primary = ivw_primary(bx, by, byse)
        for value, target in [(primary['beta'], ref.IVW_beta), (primary['se'], ref.IVW_se),
                              (primary['p'], ref.IVW_p_recalculated), (primary['q'], ref.IVW_Q)]:
            assert abs(float(value) - float(target)) <= TOL
        seed = BASE_SEED + 2000 + index
        presso = simulate_presso(bx, by, bxse, byse, seed)
        base = ref.to_dict()
        base.pop('PRESSO_global_p', None)
        base.pop('PRESSO_n_outliers', None)
        base.update({'input_sha256': sha256_bytes(raw_bytes),
                     'MR_PRESSO_implementation': 'independent Python numerical implementation, not the R package',
                     'MR_PRESSO_NbDistribution': NB, 'MR_PRESSO_seed': seed,
                     'MR_PRESSO_global_RSS_obs': presso['global_RSS_obs'],
                     'MR_PRESSO_global_p': presso['global_p'],
                     'MR_PRESSO_global_p_report': presso['global_p_report'],
                     'MR_PRESSO_outlier_test_conducted': presso['outlier_test_conducted'],
                     'MR_PRESSO_n_outliers_Bonferroni_0_05': len(presso['outlier_idx']),
                     'MR_PRESSO_outlier_SNPs': ';'.join(data.SNP.iloc[presso['outlier_idx']].astype(str)),
                     'MR_PRESSO_raw_beta': presso['raw']['beta'],
                     'MR_PRESSO_raw_se': presso['raw']['se'],
                     'MR_PRESSO_raw_p': presso['raw']['p'],
                     'MR_PRESSO_corrected_beta': presso['corrected']['beta'],
                     'MR_PRESSO_corrected_se': presso['corrected']['se'],
                     'MR_PRESSO_corrected_p': presso['corrected']['p'],
                     'MR_PRESSO_distortion_coefficient_percent': presso['distortion_coefficient_percent'],
                     'MR_PRESSO_distortion_p': presso['distortion_p']})
        reverse_rows.append(base)
        if presso['outlier_test_conducted']:
            for i, snp in enumerate(data.SNP):
                adjusted = presso['outlier_adjusted_p'][i]
                outlier_rows.append({'disease': 'winkler_earlyamd', 'layer': layer, 'SNP': snp,
                                     'RSSobs': presso['dif_obs_squared'][i],
                                     'Pvalue_raw': presso['outlier_raw_p'][i],
                                     'Pvalue_Bonferroni': adjusted,
                                     'Pvalue_report': f'<{len(data) / NB:g}' if adjusted == 0 else str(adjusted),
                                     'outlier_Bonferroni_0_05': bool(adjusted <= ALPHA)})
        for i, snp in enumerate(data.SNP):
            keep = np.arange(len(data)) != i
            fit = ivw_primary(bx[keep], by[keep], byse[keep])
            reverse_loo_rows.append({'disease': 'winkler_earlyamd', 'layer': layer,
                                     'omitted_SNP': snp, 'nsnp_remaining': int(keep.sum()),
                                     'ivw_beta': fit['beta'], 'ivw_se': fit['se'], 'ivw_p_t': fit['p'],
                                     'Q': fit['q'], 'Q_df': fit['df'], 'Q_p': fit['q_p'],
                                     'beta_change_from_full': fit['beta'] - primary['beta'],
                                     'absolute_beta_change_from_full': abs(fit['beta'] - primary['beta']),
                                     'direction_matches_full': bool(np.sign(fit['beta']) == np.sign(primary['beta']))})
            ratio = by[i] / bx[i]
            ratio_se = abs(byse[i] / bx[i])
            z = ratio / ratio_se
            wald_rows.append({'disease': 'winkler_earlyamd', 'layer': layer, 'SNP': snp,
                              'CHR': int(data.CHR.iloc[i]), 'BP': int(data.BP.iloc[i]),
                              'A1': data.A1.iloc[i], 'A2': data.A2.iloc[i],
                              'beta_exposure': bx[i], 'se_exposure': bxse[i],
                              'beta_outcome': by[i], 'se_outcome': byse[i],
                              'wald_ratio_beta': ratio, 'wald_ratio_se': ratio_se,
                              'wald_ratio_z': z, 'wald_ratio_p': 2 * stats.norm.sf(abs(z)),
                              'ci95_low': ratio - stats.norm.ppf(0.975) * ratio_se,
                              'ci95_high': ratio + stats.norm.ppf(0.975) * ratio_se,
                              'palindromic': bool(data.PALINDROMIC.iloc[i])})
reverse = pd.DataFrame(reverse_rows)
save_csv(reverse, 'reverseMR_earlyAMD_unified.csv')
save_csv(pd.DataFrame(outlier_rows), 'reverseMR_earlyAMD_MR_PRESSO_outlier_tests.csv')
save_csv(pd.DataFrame(reverse_loo_rows), 'reverseMR_earlyAMD_LOO.csv')
save_csv(pd.DataFrame(wald_rows), 'reverseMR_earlyAMD_single_SNP_Wald.csv')
assert len(reverse) == 9 and len(reverse_loo_rows) == 72 and len(wald_rows) == 72
ps = reverse[reverse.layer.eq('PS')].iloc[0]
assert ps.MR_PRESSO_outlier_test_conducted == False
assert ps.MR_PRESSO_n_outliers_Bonferroni_0_05 == 0
manifest['analyses']['reverseMR_earlyAMD'] = {
    'models': 9, 'input_IVW_exactly_reproduced_within_1e_6': True,
    'MR_PRESSO_seed_rule': '20261025 + 2000 + lexicographic layer index',
    'MR_PRESSO_NbDistribution': NB,
    'global_p': dict(zip(reverse.layer, reverse.MR_PRESSO_global_p)),
    'PS_global_p': float(ps.MR_PRESSO_global_p),
    'PS_formal_outlier_test_conducted': bool(ps.MR_PRESSO_outlier_test_conducted),
    'PS_formal_outliers': int(ps.MR_PRESSO_n_outliers_Bonferroni_0_05),
    'note': 'Historical PRESSO columns were removed from the unified table and replaced by one formal analysis.'}

manifest['validation'] = {
    'all_assertions_passed': True, 'numeric_tolerance': TOL,
    'no_uploaded_file_modified': True, 'no_external_download': True,
    'zero_simulation_p_reporting': 'Report P<0.0001 when the empirical count is zero'}
manifest_path = OUT / 'consolidated_analysis_manifest.json'
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'status': 'complete', 'output': str(OUT),
                  'analyses': list(manifest['analyses']),
                  'files': len(manifest['files'])}, ensure_ascii=False, indent=2))
