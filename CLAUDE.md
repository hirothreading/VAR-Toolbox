# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

VAR-Toolbox is a MATLAB toolbox for Vector Autoregressive (VAR) econometric analysis, authored by Ambrogio Cesa-Bianchi. It provides OLS-based VAR estimation, multiple structural identification schemes, and publication-quality plotting utilities.

## Setup

No build step required. Add the toolbox to the MATLAB path:

```matlab
addpath(genpath('/path/to/VAR-Toolbox/v3dot0/'))
```

Optional: Install [Ghostscript](https://www.ghostscript.com/) for high-quality figure export via `SaveFigure.m`.

## Repository Structure

Active code lives in `v3dot0/`. Legacy versions are in `OldVersions/` (rarely relevant).

```
v3dot0/
├── VAR/        Core estimation and structural analysis functions
├── Stats/      Summary statistics, moving averages, correlations
├── Utils/      Data manipulation and matrix utilities
├── Figure/     Plotting, date axis handling, figure export
├── Auxiliary/  Third-party algorithms (with attribution)
├── ExportFig/  Yair Altman's export_fig toolbox
├── Primer/     Tutorial script (VARToolbox_Primer.m) and sample data
└── Replic/     Replication scripts for published papers
```

## Core Architecture

### VAR Workflow

The typical workflow chains these functions:

1. **`VARmodel.m`** — Estimates a reduced-form VAR via OLS; returns a `VAR` struct containing coefficient matrices, residuals, and model metadata.
2. **`VARir.m`** / **`VARvd.m`** / **`VARhd.m`** — Compute impulse responses (IR), forecast error variance decompositions (FEVD), and historical decompositions (HD) from the `VAR` struct.
3. **`VARdrawpost.m`** — Draws bootstrap posterior distributions for IRs/FEVDs.

### Identification Methods (passed as options to `VARir.m`)

| Method | Key file(s) |
|---|---|
| Cholesky (zero short-run restrictions) | Built into `VARir.m` |
| Zero long-run restrictions | Built into `VARir.m` |
| Sign restrictions | `SR.m`, `SignRestrictions.m` |
| External instruments (proxy SVAR) | Options in `VARir.m` |
| Instruments + sign restrictions | Combined options |

### Figure Utilities

`Figure/` functions are standalone helpers used throughout:
- **`PlotSwathe.m`** — Plots confidence bands (line + shaded area)
- **`DatesCreate.m`** / **`DatesPlot.m`** — Handle date vectors for x-axes
- **`SaveFigure.m`** — Exports figures (uses ExportFig/Ghostscript when available)
- **`FigSize.m`** — Sets figure dimensions for publication

### Primer

`v3dot0/Primer/VARToolbox_Primer.m` is the canonical usage reference for linear VARs. `v3dot0/Primer/TVARToolbox_Primer.m` is the equivalent for TVARs.

---

## TVAR Extension

Eight new files in `v3dot0/VAR/` implement Threshold VARs (TVARs) and Proxy Threshold SVARs. The design principle is that `TVAR.regime{k}` is shaped identically to the `VAR` struct from `VARmodel`, so `VARir` can be delegated to for each regime without reimplementing identification logic.

### Files

| File | Signature | Purpose |
|---|---|---|
| `TVARoption.m` | `TVARopt = TVARoption` | Default options (extends `VARoption`) |
| `TVARmodel.m` | `[TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, THRVAR_EX, EXOG, nlag_ex)` | Estimate 2-regime TVAR via grid search |
| `TVARprint.m` | `TVARprint(TVAR, TVARopt, approx)` | Print threshold, regime sizes, coefficients |
| `TVARir.m` | `[IR, TVAR] = TVARir(TVAR, VARopt)` | Regime-specific IRFs; `IR` is `(nsteps × nvar × nvar × 2)` |
| `TVARirband.m` | `[INF,SUP,MED,BAR] = TVARirband(TVAR, VARopt)` | Bootstrap bands; same 4-D shape as `IR` |
| `TVARirplot.m` | `TVARirplot(IR, TVAR, VARopt, INF, SUP)` | Plot both regimes overlaid per panel |
| `TVARtest.m` | `out = TVARtest(ENDO, nlag, const, thrvar_idx, delay, VARopt, ...)` | Hansen (1999) bootstrap linearity test |
| `TVARToolbox_Primer.m` | — | End-to-end tutorial (`Primer/`) |

### TVAR Workflow

```matlab
% 1. Estimate
[TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay);
TVARopt.vnames = {'GDP','Rate'};
TVARopt.rnames = {'Recession','Expansion'};  % legend labels

% 2. Print results
TVARprint(TVAR, TVARopt);

% 3. (Optional) test linearity — Hansen (1999) bootstrap F-test
VARopt_test = TVARopt; VARopt_test.ndraws = 499; VARopt_test.method = 'wild';
out = TVARtest(ENDO, nlag, const, thrvar_idx, delay, VARopt_test);
% out.pval, out.F_stat, out.cv ([cv90 cv95 cv99]), out.grid_RSS for plotting

% 4. Regime-specific IRFs (all four identification schemes supported)
VARopt.ident = 'short';   % or 'long', 'sign', 'iv'
[IR, TVAR]   = TVARir(TVAR, VARopt);   % IR is (nsteps x nvar x nvar x 2)

% 5. Bootstrap bands
[INF, SUP, MED, BAR] = TVARirband(TVAR, VARopt);

% 6. Plot (both regimes overlaid on same axes, distinguished by colour)
TVARirplot(IR, TVAR, VARopt, INF, SUP);
```

### Proxy TVAR (`ident = 'iv'`)

Set `TVAR.IV` to a (nobs × 1) instrument vector aligned with `TVAR.ENDO` (same row count as the original ENDO before lag-trimming), then call `TVARir` as normal. The instrument is automatically trimmed to regime-specific rows using the stored `obs_idx` in each regime sub-struct.

```matlab
TVAR.IV     = my_instrument;   % (nobs x 1), NaN where unavailable
VARopt.ident = 'iv';
[IR, TVAR]  = TVARir(TVAR, VARopt);
```

### `TVAR` struct layout

`TVARmodel` returns a `TVAR` struct. Key top-level fields:

| Field | Description |
|---|---|
| `TVAR.thresh` | Estimated threshold value γ |
| `TVAR.delay` | Delay parameter d |
| `TVAR.thrvar` | `(nobse × 1)` threshold variable values q_{t-d} |
| `TVAR.regime_idx` | `(nobse × 1)` regime assignments (1 or 2) |
| `TVAR.grid_gamma` / `TVAR.grid_RSS` | Grid and RSS values for diagnostics/plotting |
| `TVAR.regime{k}` | Per-regime sub-struct (VARir-compatible) |
| `TVAR.IV` | Instrument; set by user before calling `TVARir` with `'iv'` |

Each `TVAR.regime{k}` mirrors the `VARmodel` output: `Ft`, `F`, `sigma`, `resid`, `X`, `Y`, `Fcomp`, `maxEig`, `nobs`, `B`, `Biv`, `PSI`, plus `obs_idx` (row indices into the full Y/X for IV alignment).

### `TVARoption` additional fields

Beyond the inherited `VARoption` fields:

| Field | Default | Description |
|---|---|---|
| `thrvar_idx` | `1` | Column in ENDO used as threshold variable (0 = external) |
| `delay` | `1` | Delay d for q_{t-d} |
| `trim` | `0.15` | Fraction trimmed from each tail of threshold grid (Hansen 1999) |
| `ngrid` | `300` | Grid points searched |
| `rnames` | `{'Regime 1','Regime 2'}` | Regime labels for plots and print |

### Key design decisions

- **Grid search**: searches `[quantile(q, trim), quantile(q, 1−trim)]` to ensure each regime has at least `ntotcoeff + 1` observations. Grid points that violate this constraint are skipped.
- **`TVARir` delegation**: for `'short'`/`'long'`/`'sign'`, delegates directly to `VARir(TVAR.regime{k}, VARopt)`. For `'iv'`, the ~35-line proxy SVAR identification block is run within `TVARir` using `IV_k = TVAR.IV(nlag + obs_idx, :)`, then delegates to VARir with `ident='sign'` (which trusts the pre-set `B`).
- **`'sign'` identification**: requires calling `SR(TVAR.regime{k}, SIGN, VARopt)` for each regime before `TVARir` to populate `TVAR.regime{k}.B`, matching the base toolbox convention.
- **Bootstrap (`TVARirband`)**: fixed-design — regime assignments are held at observed values when generating artificial data (Ft_k selected by original `regime_idx`); the threshold γ is re-estimated in each draw for correct coverage.
- **`TVARtest`**: bootstrap p-value computed by generating data under H₀ (linear VAR) and re-running the full grid search each draw. `out.F_boot` contains the full bootstrap distribution for plotting.
- **`TVARirplot`**: each panel shows both regimes overlaid in `cmap(1)` / `cmap(2)` with shaded confidence bands. Legend only appears on the first subplot to avoid clutter. `VARopt.rnames` (or `TVARopt.rnames`) controls the legend labels.
