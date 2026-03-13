# Change log

## TVAR Integration into VAR Toolbox v3dot0

### Overview
Added full Threshold VAR (TVAR) support to the v3dot0 VAR Toolbox. Supports both frequentist (Hansen 1997 grid search) and Bayesian (Gibbs sampler, adapted from Haroon Mumtaz) estimation, with all four identification schemes: Cholesky, long-run, sign restrictions, and external instruments (proxy SVAR).

### New files created

**Core (`v3dot0/VAR/`)**
| File | Purpose |
|------|---------|
| `TVARoption.m` | Default options (extends VARoption with TVAR-specific fields) |
| `TVARmodel.m` | TVAR estimation: freq (grid search + OLS) and bayes (Gibbs) |
| `TVARir.m` | Regime-specific IRFs (delegates to VARir per regime) |
| `TVARirband.m` | Bootstrap (freq) / posterior percentile (bayes) bands |
| `TVARirplot.m` | IRF plots: overlay (both regimes) or separate mode |
| `TVARvd.m` | Variance decomposition per regime |
| `TVARvdband.m` | VD confidence/credible bands |
| `TVARvdplot.m` | Stacked area VD plots per regime |
| `TVARhd.m` | Historical decomposition with regime-switching companion |
| `TVARhdplot.m` | HD plots with optional regime shading |
| `TVARprint.m` | Display estimation results (coefficients, eigenvalues, MCMC diag) |
| `TVARtest.m` | Hansen (1999) linearity test (updated API to use TVARopt) |

**Auxiliary (`v3dot0/Auxiliary/`)**
| File | Purpose |
|------|---------|
| `iwpQ.m` | Inverse-Wishart draw |
| `getcoef_tvar.m` | Draw VAR coefficients with stability rejection |
| `stability_tvar.m` | Companion matrix eigenvalue check |
| `loglik_tvar.m` | Gaussian log-likelihood |
| `getvarpost_tvar.m` | Log-posterior for threshold parameter |
| `create_dummies_tvar.m` | Minnesota prior dummy observations |
| `discretesample.m` | Sample from discrete distribution |
| `safeexp.m` | Numerically stable exponentiation |

**Scripts**
| File | Purpose |
|------|---------|
| `v3dot0/Replic/GK2015_TVAR/GO_GK2015_TVAR.m` | GK2015 TVAR replication |
| `v3dot0/Primer/TVARToolbox_Primer.m` | TVAR tutorial script |

### Key design decisions

1. **`TVAR.regime{k}` mirrors VAR struct**: Each regime struct has the same fields as VARmodel output (Ft, F, sigma, Fcomp, resid, Y, X, etc.), so existing functions (VARir, VARvd) work directly on each regime.

2. **Fixed-regime bootstrap**: For frequentist bands, residuals are resampled within each regime with the threshold/regime assignment held fixed. This prevents the "garbled IRFs" problem from the previous attempt.

3. **Adaptive stability threshold**: Bootstrap stability threshold is set to `max(1.05, maxEig_pt + 0.10)` where `maxEig_pt` is the point estimate's maximum eigenvalue. This prevents systematic bias from rejecting draws near the point estimate.

4. **IV pairs bootstrap**: For proxy SVAR identification, the bootstrap uses residual-IV pairs resampling (same indices for both), since the wild bootstrap cancels when applied simultaneously to residuals and instruments with fixed coefficients.

5. **Coefficient ordering**: Mumtaz's Gibbs sampler uses [lags|constant] ordering internally; converted to VARmodel's [constant|lags] ordering when building regime structs.

### Test results (GK2015 data, Cholesky, 200 draws)
- Regime 1 (308 obs, Low EBP): IR within 68% band at 96% of steps
- Regime 2 (76 obs, High EBP): IR within 68% band at 88% of steps
- VD rows sum to 100%
- HD reconstruction error: 1.78e-15
- All plotting functions work correctly

### Lessons / patterns
- **Stability threshold must match point estimate**: A hardcoded threshold below the point estimate's maxEig creates systematic downward bias (all draws near the point estimate get rejected).
- **Wild bootstrap doesn't work for IV with fixed-X**: Sign-flipping residuals and IV simultaneously with rr^2=1 cancels in the first-stage regression. Must use pairs bootstrap instead.
- **Jensen's inequality**: Bootstrap mean of impact response (sqrt of sigma) is systematically below the point estimate. This is a known theoretical property, not a bug.
- **Complex IRFs in small regimes**: IV identification can produce complex values when regime has too few observations relative to parameters. Added warning in TVARir.

---

## TVAR IRF Fix: Linear IRFs → Generalized IRFs (KPP 1996)

### Problem
The original `TVARir.m` delegated to `VARir()` per regime, computing **regime-specific linear IRFs** using the Wold representation (companion matrix). This is incorrect for a nonlinear TVAR model because:
1. The regime is held fixed throughout the IRF horizon — no regime switching possible
2. Historical starting conditions are ignored
3. The resulting IRFs do not capture the nonlinear dynamics that motivate using a TVAR

### Root cause analysis
Comparison with the reference Mumtaz code (`TVAR/getimpulse.m`) revealed it computes **Generalized IRFs (GIRFs)** via forward simulation. However, the Mumtaz code itself contains a bug when `transform=0`: it evaluates `yhatm(fi,:) = yhat(fi,:)` **before** `yhat(fi,:)` is computed, making `ystar(fi) = 0` always. Since the threshold is typically positive (tar ≈ 0.93), regime 1 is always selected during simulation, effectively reducing the GIRF to a regime-1 linear IRF. The alternative function `get_irfTVAR.m` correctly uses `yhat(fi-delay, tvar)` but is not called in the Mumtaz workflow.

### Solution
Implemented proper GIRFs following Koop, Pesaran & Potter (1996):

**New file: `TVARgirf.m`**
- For each regime, for each historical starting point in that regime:
  - Simulates unshocked and shocked paths forward, allowing endogenous regime switching
  - Uses the SAME random innovations for both paths (variance reduction)
  - Correctly evaluates threshold as `yhat(fi-delay, thrvar_idx)` (the lagged simulated value)
  - GIRF = average of (shocked - unshocked) across simulations and starting points

**Updated files:**
| File | Change |
|------|--------|
| `TVARgirf.m` | **NEW**: Core GIRF computation |
| `TVARir.m` | Rewritten to call TVARgirf instead of VARir per regime |
| `TVARirband.m` | Bayesian: GIRF at each posterior draw. Freq: GIRF at each bootstrap draw |
| `TVARoption.m` | Added `nreps_girf` (sims per starting point) and `shock_scale` |
| `GO_GK2015_TVAR.m` | Switched to Bayesian estimation + GIRF workflow |
| `TVARToolbox_Primer.m` | Updated to match new API |

### Improvements over Mumtaz reference
1. **Correct threshold evaluation**: Uses `yhat(fi-delay, thrvar_idx)` — the properly lagged simulated value — enabling endogenous regime switching during the IRF horizon
2. **Variance reduction**: Same random innovations for shocked and unshocked paths, so noise cancels in the difference
3. **Integrated into Toolbox API**: Works with TVARmodel output, posterior draws, and existing plotting infrastructure

### Test results (Mumtaz data, 4 vars, 4 lags, Bayesian)
- GIRFs converge back to zero (unlike the previous linear IRFs which could diverge)
- Economically sensible responses: credit shock → GDP contraction, T-bill decline, spread returns to zero
- Regime-specific differences: Regime 2 (high spread) shows larger initial impact
- Bayesian credible bands have correct coverage

### Test results (GK2015 data, 4 vars, 12 lags, Bayesian)
- Threshold ≈ 0.38 (EBP level), R1: 299 obs (Low EBP), R2: 85 obs (High EBP)
- Acceptance rate ≈ 0.20 (reasonable for MH)
- Regime 2 EBP impact (0.32) > Regime 1 (0.20) — financial stress amplifies credit shocks
- GDP decline larger in Regime 2 (-0.008) vs Regime 1 (-0.003)
- All plotting and band functions work correctly

### Current limitations
- GIRF identification limited to Cholesky ('short') — sign restrictions and IV not yet supported for GIRFs
- GIRFs are computationally intensive (forward simulation × starting points × posterior draws)
- VD and HD still use linear regime-specific approach (no GIRF analogue)

### Lessons
- **Always verify against reference code line-by-line**: The conceptual difference (linear vs GIRF) was the root cause, but only detailed comparison revealed the Mumtaz bug
- **Bugs in reference code can mask problems**: The Mumtaz threshold evaluation bug made its "GIRFs" equivalent to linear IRFs, which happened to look reasonable. The correct GIRF produces qualitatively similar but quantitatively different results
- **Variance reduction matters**: Using the same random draws for shocked/unshocked paths is essential for clean GIRFs with reasonable number of simulations
