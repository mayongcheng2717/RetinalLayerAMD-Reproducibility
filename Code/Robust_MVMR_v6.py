#!/usr/bin/env python3
"""Robust early-AMD MVMR sensitivity using published MVMR-Median and MVMR-Egger.

No downloads and no modification of uploaded files. MVMR-Median follows Grant and
Burgess (Statistics in Medicine 2021, DOI 10.1002/sim.9156): weighted L1 regression
with outcome-SE weights and parametric bootstrap including exposure and outcome
association uncertainty. The linear program is an independent SciPy implementation.
"""
from pathlib import Path
import hashlib, json, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.optimize import linprog
import statsmodels.api as sm

ROOT = Path('/mnt/agents/output/审稿修订_本地补分析包')
CON = ROOT / 'work_v6' / 'consolidated_results'
MATRIX = ROOT / 'returned_results/02_mvmr_ucsc_verified/matrices/winkler_earlyamd.tsv.gz'
LAYERS = ['INL', 'OPL', 'RPE']
BOOT = 10000
SEED = 20261026

mm = pd.read_csv(MATRIX, sep='\t')
BX = mm[[f'bx_{x}' for x in LAYERS]].to_numpy(float)
BXSE = mm[[f'bxse_{x}' for x in LAYERS]].to_numpy(float)
by = mm.by.to_numpy(float)
byse = mm.byse.to_numpy(float)


def weighted_median_fit(X, y, outcome_se):
    """Minimize sum outcome_se^-2 * abs(y-X beta), no intercept."""
    n, k = X.shape
    weights = 1 / outcome_se**2
    # y = X beta + u - v, with u,v >= 0 and beta unconstrained.
    objective = np.r_[np.zeros(k), weights, weights]
    equality = np.c_[X, np.eye(n), -np.eye(n)]
    bounds = [(None, None)] * k + [(0, None)] * (2*n)
    result = linprog(objective, A_eq=equality, b_eq=y, bounds=bounds, method='highs')
    if not result.success:
        raise RuntimeError(result.message)
    return result.x[:k], float(result.fun)


def egger_fit(X, y, outcome_se):
    # Orientate on the first exposure, as in mr_mvegger(orientate=1).
    orientation = np.sign(X[:, 0])
    Xo, yo = X * orientation[:, None], y * orientation
    design = np.c_[np.ones(len(yo)), Xo]
    weights = 1 / outcome_se**2
    sw = np.sqrt(weights)
    theta = np.linalg.lstsq(design * sw[:, None], yo * sw, rcond=None)[0]
    residual = yo - design @ theta
    q = float(np.sum(weights * residual**2))
    df = len(yo) - design.shape[1]
    phi = max(q / df, 1.0)
    covariance = np.linalg.pinv(design.T @ (weights[:, None] * design)) * phi
    se = np.sqrt(np.diag(covariance))
    p = 2 * stats.t.sf(np.abs(theta / se), df)
    critical = stats.t.ppf(0.975, df)
    return theta, se, p, theta-critical*se, theta+critical*se, q, df, phi

# Verify the LP against an independent unweighted quantile-regression implementation.
lp_beta, lp_objective_unweighted = weighted_median_fit(BX, by, np.ones(len(by)))
qr_fit = sm.QuantReg(by, BX).fit(q=0.5)
qr_beta = qr_fit.params
qr_objective_unweighted = float(np.sum(np.abs(by-BX@qr_beta)))
lp_check = float(np.max(np.abs(lp_beta-qr_beta)))
objective_check = float(lp_objective_unweighted-qr_objective_unweighted)
# QuantReg uses an iterative numerical solver. The LP is accepted when estimates
# agree to 1e-4 and its L1 objective is no worse within numerical tolerance.
if lp_check > 1e-4 or objective_check > 1e-6:
    raise RuntimeError(f'Independent quantile regression check failed: beta={lp_check}, objective={objective_check}')

rng = np.random.default_rng(SEED)
summary_rows = []
bootstrap_frames = []
models = [('all_118_instruments', np.arange(len(mm))),
          ('without_rs10737680', np.flatnonzero(~mm.SNP.eq('rs10737680').to_numpy()))]
for model_name, indices in models:
    X, XSE = BX[indices], BXSE[indices]
    y, yse = by[indices], byse[indices]
    point, objective = weighted_median_fit(X, y, yse)
    boot = np.empty((BOOT, len(LAYERS)))
    for iteration in range(BOOT):
        xboot = rng.normal(X, XSE)
        yboot = rng.normal(y, yse)
        boot[iteration] = weighted_median_fit(xboot, yboot, yse)[0]
    boot_se = boot.std(axis=0, ddof=1)
    for j, layer in enumerate(LAYERS):
        z = point[j] / boot_se[j]
        summary_rows.append({
            'outcome': 'winkler_earlyamd', 'model': model_name, 'method': 'MVMR-Median',
            'layer': layer, 'nsnp': len(indices), 'beta': point[j], 'se_parametric_bootstrap': boot_se[j],
            'ci95_low': point[j]-stats.norm.ppf(0.975)*boot_se[j],
            'ci95_high': point[j]+stats.norm.ppf(0.975)*boot_se[j],
            'p_normal': 2*stats.norm.sf(abs(z)),
            'bootstrap_percentile_ci95_low': np.quantile(boot[:, j], 0.025),
            'bootstrap_percentile_ci95_high': np.quantile(boot[:, j], 0.975),
            'bootstrap_same_sign_fraction': float(np.mean(np.sign(boot[:, j]) == np.sign(point[j]))),
            'bootstrap_iterations': BOOT, 'seed': SEED,
            'weighted_L1_objective': objective,
            'interpretation': 'Pleiotropy-robust sensitivity; does not model sample overlap or same-sample OCT exposure covariance'})
        bf = pd.DataFrame({'model': model_name, 'iteration': np.arange(1, BOOT+1),
                           'layer': layer, 'beta_bootstrap': boot[:, j]})
        bootstrap_frames.append(bf)

median = pd.DataFrame(summary_rows)
median.to_csv(CON / 'MVMR_earlyAMD_median_sensitivity.csv', index=False)
pd.concat(bootstrap_frames, ignore_index=True).to_csv(CON / 'MVMR_earlyAMD_median_bootstrap.csv.gz', index=False, compression='gzip')

# MVMR-Egger sensitivity on all 118 instruments, with coherent t inference df=n-k-1.
theta, se, p, low, high, q, df, phi = egger_fit(BX, by, byse)
egger_rows = []
for i, term in enumerate(['Intercept'] + LAYERS):
    egger_rows.append({'outcome': 'winkler_earlyamd', 'method': 'MVMR-Egger', 'term': term,
                       'nsnp': len(mm), 'beta': theta[i], 'se': se[i], 'df': df,
                       'p_t': p[i], 'ci95_low': low[i], 'ci95_high': high[i],
                       'Q': q, 'Q_df': df, 'Q_p': stats.chi2.sf(q, df),
                       'residual_scale_inflation': phi,
                       'orientation': 'sign of INL exposure association; orientate=1',
                       'interpretation': 'Directional-pleiotropy sensitivity under InSIDE; not a resolution of sample overlap or model misspecification'})
egger = pd.DataFrame(egger_rows)
egger.to_csv(CON / 'MVMR_earlyAMD_Egger_sensitivity.csv', index=False)

manifest = {
    'status': 'complete', 'method_reference': 'Grant AJ, Burgess S. Stat Med. 2021;40:5813-5830. DOI 10.1002/sim.9156',
    'method_source': 'https://github.com/aj-grant/robust-mvmr/blob/master/mvmr_functions.R',
    'source_file_sha_verified_from_GitHub_API': '4e3b1ea32ba9da971e24ca248bb73cdf5233423f',
    'implementation': 'Independent SciPy HiGHS weighted-L1 linear program and NumPy weighted least squares; no R package execution claimed',
    'bootstrap_iterations': BOOT, 'seed': SEED,
    'bootstrap_draws': 'Independent normal draws for every SNP-exposure and SNP-outcome association using reported SE; fixed original outcome-SE weights',
    'inference': 'MVMR-Median bootstrap SE with normal P and CI; MVMR-Egger t P and matched t CI with df=n-(k+1)',
    'LP_crosscheck_max_abs_beta_difference_vs_statsmodels_unweighted': lp_check,
    'software': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
    'limitations': ['Actual participant overlap is unknown and not modelled',
                    'Same-sample OCT exposure covariance is not modelled',
                    'Robust estimators do not repair broad model misspecification or prove causality'],
    'files': {}}
for name in ['MVMR_earlyAMD_median_sensitivity.csv', 'MVMR_earlyAMD_median_bootstrap.csv.gz',
             'MVMR_earlyAMD_Egger_sensitivity.csv']:
    path = CON / name
    manifest['files'][name] = {'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
(CON / 'MVMR_earlyAMD_robust_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
print(median[['model', 'layer', 'beta', 'se_parametric_bootstrap', 'p_normal', 'ci95_low', 'ci95_high']].to_string(index=False))
print(egger[['term', 'beta', 'se', 'p_t', 'Q', 'Q_df', 'Q_p']].to_string(index=False))
print(json.dumps({'status': 'complete', 'LP_crosscheck': lp_check}, indent=2))
