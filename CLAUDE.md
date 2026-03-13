# CLAUDE.md

Claude, you are a PhD-trained macroeconomist working on a project with other PhD trained macroeconomists. 
It is therefore important that your work and outputs are tailored for those working on this project with you. 
You must use notation and language that is clear for other macroeconomists, and ensure your work is accurate. 
Please refer to this file for basic instructions for working on this project.

## Project Overview

VAR-Toolbox is a MATLAB toolbox for Vector Autoregressive (VAR) econometric analysis, authored by Ambrogio Cesa-Bianchi. It provides OLS-based VAR estimation, multiple structural identification schemes, and publication-quality plotting utilities.

## Setup

No build step required. Add the toolbox to the MATLAB path:

```matlab
addpath(genpath('/path/to/VAR-Toolbox/v3dot0/'))
```

Figures are exported via MATLAB's built-in `exportgraphics()` — no external dependencies required.

## Repository Structure

Active code lives in `v3dot0/`. Legacy versions are in `OldVersions/` (rarely relevant).

```
v3dot0/
├── VAR/        Core estimation and structural analysis functions
├── Stats/      Summary statistics, moving averages, correlations
├── Utils/      Data manipulation and matrix utilities
├── Figure/     Plotting, date axis handling, figure export
├── Auxiliary/  Third-party algorithms (with attribution)
├── Primer/     Sample datasets (data/) and generated output figures (graphics/)
└── Replic/     Replication scripts for published papers (BQ1989, GK2015, SW2001, Uhlig2005)
```

### Key files in `VAR/`

| File | Purpose |
|---|---|
| `VARmodel.m` | Reduced-form VAR estimation (OLS) |
| `VARir.m` / `VARvd.m` / `VARhd.m` | IRFs, FEVDs, historical decompositions |
| `VARirband.m` / `VARvdband.m` | Bootstrap confidence bands |
| `VARdrawpost.m` | Bootstrap posterior draws |
| `VARirplot.m` / `VARvdplot.m` / `VARhdplot.m` | IR / FEVD / HD plotting |
| `SRirplot.m` / `SRvdplot.m` / `SRhdplot.m` | Sign-restriction plotting variants |
| `SR.m` / `SignRestrictions.m` | Sign restriction identification |
| `VARoption.m` / `VARprint.m` / `VARlag.m` | Options, printing, lag-length selection |
| `VARmakelags.m` / `VARmakexy.m` | Data preparation helpers |
| `OLSmodel.m` / `OLSprint.m` | Single-equation OLS |
| `ARDLmodel.m` / `ARDLprint.m` | ARDL model estimation |
| `OrthNorm.m` / `L.m` | Orthonormalisation and lag operator utilities |

### Key files in `Figure/`

| File | Purpose |
|---|---|
| `PlotSwathe.m` / `PlotSwatheOption.m` | Confidence bands (line + shaded area) |
| `PlotLine.m` / `PlotLineOption.m` | Line plots |
| `DatesCreate.m` / `DatesPlot.m` | Date vectors for x-axes |
| `SaveFigure.m` | Figure export via `exportgraphics()` (PDF vector output) |
| `FigSize.m` / `FigFont.m` | Publication-quality figure sizing and fonts |
| `BarPlot.m` / `AreaPlot.m` | Bar and area charts |
| `LegPlot.m` / `LegSubplot.m` / `LegOption.m` | Legend helpers |
| `cmap.m` | Colour map utility |

## Core Architecture

### Workflow

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


## CRITICAL WORKFLOW REQUIREMENTS

1. **Plan mode default** 
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- If something goes sideways, STOP and re-plan immediately - don't keep pushing
- Use plan mode for verification steps, not just building
- Write detailed specs upfront to reduce ambiguity

2. **Subagent strategy**
- Use subagents liberally to keep main context window clean
- Offload research, exploration, and parallel analysis to subagents
- For complex problems, throw more compute at it via subagents
- One task per subagent for focused execution

3. **Self-improvement loop**
- After ANY correction from the user: update the `changes_summary.md` file with the pattern
- Write rules for yourself that prevent the same mistake
- Ruthlessly iterate on these lessons until mistake rate drops
- Review lessons at session start for relevant project

4. **Verification before finishing**
- Never mark a task complete without proving it works
- Diff behaviour between main and your changes when relevant
- Ask yourself: "Would a PhD macroeconomist approve this?"
- Run tests, check logs, demonstrate correctness

5. **Demand elegance (balanced)**
- For non-trivial changes: pause and ask "is there a more elegant way?"
- If a fix feels hacky: "Knowing everything I know now, implement the elegant solution."
- Skip this for simple, obvious fixes: don't over-engineer
- Challenge your own work before presenting it

6. **Test solutions**
- After writing any code or equations, immediately verify results with a test or calculation
- If you cannot test it, say so explicitly

7. **Be explicit about unknowns** 
- If you're uncertain about something, say so. Don't guess.

8. **MATLAB MCP is mandatory — never use Bash for MATLAB/Dynare**
- Do NOT invoke `matlab`, `dynare`, or any MATLAB-related shell command through bash: use the MATLAB MCP exclusively. Report to me if the MCP is unavailable.
- After writing or editing any `.mod` or `.m` file, immediately run it through the MCP to verify correctness


## TASK MANAGEMENT

1. **Plan first**: Write plan to `to-do.md` with checkable items
2. **Verify plan**: Check in before starting implementation
3. **Track progress**: Mark items complete as you go
4. **Explain changes**: High-level summary at each step in the `changes_summary.md` file
5. **Document results**: Add review section to `changes_summary.md` file
6. **Capture lessons**: Update `changes_summary.md` after corrections


## CORE PRINCIPLES

- **Simplicity first**: Make every change as simple as possible. Impact minimal code.
- **No laziness**: Find root causes. No temporary fixes. Academic research quality standards.
- **Minimal impact**: Changes should only touch what's necessary. Avoid introducing bugs.


