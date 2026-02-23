%% TVAR TOOLBOX PRIMER
%==========================================================================
% This script illustrates the main features of the TVAR extension to
% VAR-Toolbox v3.0. It covers:
%
%   1. Data loading and treatment
%   2. TVAR estimation (TVARmodel)
%   3. Printing results (TVARprint)
%   4. Threshold linearity test (TVARtest)
%   5. Plotting the RSS objective function
%   6. Regime-specific IRFs with Cholesky identification (TVARir, TVARirband)
%   7. Proxy Threshold SVAR IRFs (TVARir with ident='iv')
%
% The same US macro dataset used in VARToolbox_Primer.m is used here
% (Simple_Data.xlsx: real GDP growth and 1-year Treasury Bill rate).
%
% Identification interpretation:
%   Regime 1 (Low growth)  : GDP growth <= estimated threshold
%   Regime 2 (High growth) : GDP growth >  estimated threshold
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================

clear; close all; clc
warning off all
format short g
addpath(genpath('../'))   % add v3dot0 to path if running from Primer folder


%% 1. LOAD AND TREAT DATA
%--------------------------------------------------------------------------
% Load the Simple_Data dataset: quarterly US macro data.
[xlsdata, xlstext] = xlsread('data/Simple_Data.xlsx','Sheet1');
dates      = xlstext(3:end, 1);
datesnum   = Date2Num(dates);
vnames_long = xlstext(1, 2:end);
vnames      = xlstext(2, 2:end);
nvar_raw    = length(vnames);
data        = Num2NaN(xlsdata);

for ii = 1:nvar_raw
    DATA.(vnames{ii}) = data(:,ii);
end
nobs = size(data,1);

% Compute log-differences and first differences
tempnames = {'gdp','cpi','i1yr'};
temptreat = {'logdiff','logdiff','diff'};
tempscale = [100, 100, 1];
for ii = 1:length(tempnames)
    aux = {['d' tempnames{ii}]};
    DATA.(aux{1}) = tempscale(ii) * XoX(DATA.(tempnames{ii}), 1, temptreat{ii});
end


%% 2. TVAR ESTIMATION
%--------------------------------------------------------------------------
% We use a 2-variable VAR (GDP growth and the 1-year rate) with 2 lags
% and a constant. The threshold variable is GDP growth (column 1), lagged
% by 1 period, so regime switches when gdp growth_{t-1} crosses a threshold.

Xvnames      = {'dgdp','i1yr'};
Xvnames_long = {'Real GDP Growth','1-yr Interest Rate'};
Xnvar        = length(Xvnames);

X = nan(nobs, Xnvar);
for ii = 1:Xnvar
    X(:,ii) = DATA.(Xvnames{ii});
end
[X, fo, lo] = CommonSample(X);

% VAR settings
det   = 1;    % constant
nlags = 2;    % lag order

% Estimate the TVAR
%   thrvar_idx = 1  => threshold variable is column 1 (GDP growth)
%   delay      = 1  => use q_{t-1}
[TVAR, TVARopt] = TVARmodel(X, nlags, det, 1, 1);

% Set labels
TVARopt.vnames = Xvnames_long;
TVARopt.rnames = {'Low growth (Regime 1)', 'High growth (Regime 2)'};


%% 3. PRINT RESULTS
%--------------------------------------------------------------------------
TVARprint(TVAR, TVARopt, 3);


%% 4. THRESHOLD LINEARITY TEST (Hansen 1999)
%--------------------------------------------------------------------------
% Bootstrap test of H0: linear VAR  vs.  H1: two-regime TVAR.
% Uses 499 bootstrap draws.

VARopt_test          = TVARopt;
VARopt_test.ndraws   = 499;
VARopt_test.method   = 'wild';   % wild bootstrap (more robust to heteroskedasticity)
VARopt_test.mult     = 100;

out_test = TVARtest(X, nlags, det, 1, 1, VARopt_test);

fprintf('\n  Conclusion: p-value = %.3f ', out_test.pval);
if out_test.pval < 0.05
    fprintf('=> Reject H0 at 5%% (evidence of threshold nonlinearity)\n\n');
else
    fprintf('=> Cannot reject H0 at 5%%\n\n');
end


%% 5. PLOT THE RSS OBJECTIVE FUNCTION
%--------------------------------------------------------------------------
% The grid_RSS plot shows how total RSS varies with the threshold candidate.
% The minimum identifies the estimated threshold.

FigSize(26,8);
subplot(1,2,1)
plot(out_test.grid_gamma, out_test.grid_RSS, 'LineWidth', 2, 'Color', cmap(1));
hold on;
xline(out_test.thresh, '--k', 'LineWidth', 1.5);
xlabel('Threshold candidate (\gamma)');
ylabel('Total RSS');
title('RSS objective function');
grid on;

subplot(1,2,2)
% Bootstrap distribution of F-statistic
histogram(out_test.F_boot, 30, 'FaceColor', cmap(2), 'EdgeColor','w'); hold on;
xline(out_test.F_stat, '-k', 'LineWidth', 2);
xlabel('F-statistic');
title(['Bootstrap distribution | F = ' num2str(out_test.F_stat,'%.2f') ...
       ' | p = ' num2str(out_test.pval,'%.3f')]);
grid on;

SaveFigure('graphics/TVAR_test', 2);
clf('reset');


%% 6. CHOLESKY IRFs WITH BOOTSTRAP CONFIDENCE BANDS
%--------------------------------------------------------------------------
% Identification: zero contemporaneous restrictions (Cholesky, lower-
% triangular). First shock = demand shock (GDP ordered first).

VARopt_chol                = TVARopt;
VARopt_chol.ident          = 'short';
VARopt_chol.nsteps         = 20;
VARopt_chol.ndraws         = 500;
VARopt_chol.pctg           = 68;
VARopt_chol.method         = 'bs';
VARopt_chol.mult           = 100;
VARopt_chol.FigSize        = [26, 16];
VARopt_chol.figname        = 'graphics/TVAR_chol_';
VARopt_chol.quality        = 2;
VARopt_chol.suptitle       = 0;
VARopt_chol.snames         = {'\epsilon^{GDP}', '\epsilon^{Rate}'};

% Point estimates
[IR_chol, TVAR] = TVARir(TVAR, VARopt_chol);

% Bootstrap confidence bands
[INF_chol, SUP_chol, MED_chol, BAR_chol] = TVARirband(TVAR, VARopt_chol);

% Plot: both regimes overlaid, with shaded 68% bands
TVARirplot(IR_chol, TVAR, VARopt_chol, INF_chol, SUP_chol);


%% 7. PROXY THRESHOLD SVAR (external instrument identification)
%--------------------------------------------------------------------------
% Identification: external instrument for the GDP growth shock.
% We construct an artificial instrument (demand shock from Cholesky
% + noise) as a stand-in for a real instrument.
%
% In practice, replace 'iv' with a real external instrument (e.g.,
% a measure of fiscal or monetary policy shocks).

rng(1);
% Extract Cholesky demand shock from both regimes' residuals, reassembled
eps_chol = nan(TVAR.nobs, Xnvar);
eps_chol(TVAR.regime_idx==1, :) = (TVAR.regime{1}.B \ TVAR.regime{1}.resid')';
eps_chol(TVAR.regime_idx==2, :) = (TVAR.regime{2}.B \ TVAR.regime{2}.resid')';

noise = randn(size(X,1), 1);
noise = noise(1:TVAR.nobs);            % match effective sample
iv_raw = [NaN; eps_chol(1:end-1, 1)] + 0.5*noise;   % lagged demand shock + noise

% Pad back to full ENDO length for TVARir (must match TVAR.ENDO)
iv_full = nan(size(TVAR.ENDO, 1), 1);
% TVAR.Y rows correspond to ENDO rows nlag+1:end
iv_full(nlags+1:end) = iv_raw;

TVAR.IV = iv_full;

VARopt_iv                  = TVARopt;
VARopt_iv.ident            = 'iv';
VARopt_iv.nsteps           = 20;
VARopt_iv.ndraws           = 500;
VARopt_iv.pctg             = 68;
VARopt_iv.method           = 'wild';
VARopt_iv.mult             = 100;
VARopt_iv.FigSize          = [26, 16];
VARopt_iv.figname          = 'graphics/TVAR_iv_';
VARopt_iv.quality          = 2;
VARopt_iv.suptitle         = 1;
VARopt_iv.snames           = {'\epsilon^{GDP}', '\epsilon^{Other}'};
VARopt_iv.rnames           = {'Low growth (Regime 1)', 'High growth (Regime 2)'};

% Point estimates
[IR_iv, TVAR] = TVARir(TVAR, VARopt_iv);

% Bootstrap bands
[INF_iv, SUP_iv, ~, ~] = TVARirband(TVAR, VARopt_iv);

% Plot
TVARirplot(IR_iv, TVAR, VARopt_iv, INF_iv, SUP_iv);


%% 8. SUMMARY: SANITY CHECKS
%--------------------------------------------------------------------------
fprintf('\n===== Sanity checks =====\n');

% Check regime obs sum to total
assert(TVAR.regime{1}.nobs + TVAR.regime{2}.nobs == TVAR.nobs, ...
    'Regime obs do not sum to total!');
fprintf('  Regime obs: %d + %d = %d  [OK]\n', ...
    TVAR.regime{1}.nobs, TVAR.regime{2}.nobs, TVAR.nobs);

% Cholesky: B*B'' should equal sigma within each regime
for k = 1:2
    B_k     = TVAR.regime{k}.B;
    sigma_k = TVAR.regime{k}.sigma;
    err_k   = max(max(abs(B_k*B_k' - sigma_k)));
    fprintf('  Regime %d  B*B'' - sigma max abs err: %.2e  ', k, err_k);
    if err_k < 1e-10; fprintf('[OK]\n'); else; fprintf('[CHECK]\n'); end
end

% Residuals reassembled from regimes should match TVAR.Y
resid_full = nan(TVAR.nobs, TVAR.nvar);
resid_full(TVAR.regime_idx==1,:) = TVAR.regime{1}.resid;
resid_full(TVAR.regime_idx==2,:) = TVAR.regime{2}.resid;
% Reconstruct fitted values per regime
Yhat_full = nan(TVAR.nobs, TVAR.nvar);
Yhat_full(TVAR.regime_idx==1,:) = TVAR.X(TVAR.regime_idx==1,:) * TVAR.regime{1}.Ft;
Yhat_full(TVAR.regime_idx==2,:) = TVAR.X(TVAR.regime_idx==2,:) * TVAR.regime{2}.Ft;
recon_err = max(max(abs(TVAR.Y - (Yhat_full + resid_full))));
fprintf('  Y reconstruction error: %.2e  ', recon_err);
if recon_err < 1e-10; fprintf('[OK]\n'); else; fprintf('[CHECK]\n'); end

fprintf('=========================\n\n');

close all;
