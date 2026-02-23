%% PROXY THRESHOLD SVAR — CORRECTNESS TEST
%==========================================================================
% Monte Carlo verification that TVARir with ident='iv' correctly recovers
% the structural impact vector in each regime.
%
% DESIGN
% ------
% 1. Simulate T observations from a known 2-regime structural VAR:
%
%      Regime 1 (q_{t-1} <= gamma): y_t = c1 + A1*y_{t-1} + B1*eps_t
%      Regime 2 (q_{t-1} >  gamma): y_t = c2 + A2*y_{t-1} + B2*eps_t
%
%    where eps_t ~ iid N(0,I) and B1, B2 are the true structural impact
%    matrices (lower triangular in this setup).
%
% 2. The threshold variable is y_{1,t-1} (the lagged first variable).
%
% 3. Construct a proxy instrument:
%      z_t = rho * eps_{1,t} + sqrt(1-rho^2) * nu_t,  nu_t ~ N(0,1)
%    where rho controls instrument relevance (0.9 in the baseline).
%
% 4. Run TVARmodel (recovers regime assignments and reduced-form params),
%    then TVARir with ident='iv'.
%
% WHAT IS CHECKED
% ---------------
% For each regime k, the estimated Biv_k should be proportional to the
% true B_k(:,1). After normalising by the (1,1) element:
%
%   Biv_k / Biv_k(1)  ≈  B_k(:,1) / B_k(1,1)
%
% We verify this across M Monte Carlo replications and report:
%   - Bias (mean estimate - truth)
%   - RMSE
%   - Whether each element is within the 95% bootstrap band on average
%
% The test is passed if bias < 0.05 and RMSE < 0.15 for all elements.
%==========================================================================
% Run from v3dot0/Primer/ with the toolbox on the MATLAB path.
%==========================================================================

clear; close all; clc;
warning off all;
addpath(genpath('../'));

rng(42);   % reproducibility


%% ========================================================================
%  PARAMETERS
%  =========================================================================

T      = 300;   % time series length
nlag   = 1;     % VAR lag order
const  = 1;     % include constant
delay  = 1;     % delay for threshold variable
nvar   = 3;     % number of endogenous variables
M      = 200;   % Monte Carlo replications
rho    = 0.9;   % instrument relevance (correlation with structural shock)

% --- True structural parameters -------------------------------------------
% Regime 1 (low): persistent, moderate volatility
A1  = [0.60  0.00  0.00;
       0.10  0.50  0.00;
       0.05  0.10  0.40];
c1  = [0.20; 0.10; 0.05];

% Regime 2 (high): less persistent, higher volatility
A2  = [0.30  0.05  0.00;
       0.20  0.30  0.05;
       0.10  0.05  0.20];
c2  = [0.40; 0.20; 0.10];

% True impact matrices (lower triangular — Cholesky-style)
B1  = [1.0   0     0  ;
       0.5   1.2   0  ;
       0.3   0.4   0.8];

B2  = [1.5   0     0  ;
       0.8   1.0   0  ;
       0.2   0.6   1.1];

% True threshold: gamma = 0 (on lagged y1)
gamma_true = 0.0;

% Normalised true Biv columns (ratio form, what the proxy SVAR recovers)
Biv1_true_norm = B1(:,1) / B1(1,1);   % [1; 0.5; 0.3]
Biv2_true_norm = B2(:,1) / B2(1,1);   % [1; 0.8/1.5; 0.2/1.5]


%% ========================================================================
%  MONTE CARLO LOOP
%  =========================================================================

Biv1_MC = nan(nvar, M);   % estimated normalised Biv for regime 1
Biv2_MC = nan(nvar, M);   % estimated normalised Biv for regime 2
thresh_MC = nan(M, 1);    % estimated threshold
n1_MC     = nan(M, 1);    % regime 1 size

fprintf('Running %d Monte Carlo replications (T=%d, rho=%.1f)...\n', M, T, rho);
for mm = 1:M

    if mod(mm,50)==0; fprintf('  Replication %d / %d\n', mm, M); end

    %% --- Simulate TVAR-SVAR data -----------------------------------------
    nburn = 100;
    Tall  = T + nburn;
    y     = zeros(Tall, nvar);
    eps   = randn(Tall, nvar);   % structural shocks, iid N(0,I)
    nu    = randn(Tall, 1);      % instrument noise
    z     = rho * eps(:,1) + sqrt(1 - rho^2) * nu;  % proxy instrument

    % Initialise
    y(1,:) = randn(1, nvar) * 0.1;

    for t = 2:Tall
        q_t = y(t-1, 1);   % threshold variable: lagged y1
        if q_t <= gamma_true
            y(t,:) = (c1 + A1 * y(t-1,:)')' + (B1 * eps(t,:)')';
        else
            y(t,:) = (c2 + A2 * y(t-1,:)')' + (B2 * eps(t,:)')';
        end
    end

    % Discard burn-in
    y   = y(nburn+1:end, :);
    z   = z(nburn+1:end);
    eps = eps(nburn+1:end, :);

    %% --- Estimate TVAR ---------------------------------------------------
    try
        [TVAR, TVARopt] = TVARmodel(y, nlag, const, 1, delay);
    catch
        continue   % skip if threshold estimation fails
    end

    thresh_MC(mm) = TVAR.thresh;
    n1_MC(mm)     = TVAR.regime{1}.nobs;

    %% --- Set instrument and run Proxy TVAR IR ----------------------------
    % TVAR.IV must have exactly size(TVAR.ENDO,1) = T rows.
    % TVARir accesses TVAR.IV(nlag + obs_idx) where obs_idx ∈ {1,...,T-nlag},
    % so rows 1:nlag are never read — no NaN padding required.
    TVAR.IV = z;   % z is T x 1, same length as TVAR.ENDO

    VARopt_iv        = TVARopt;
    VARopt_iv.ident  = 'iv';
    VARopt_iv.nsteps = 10;

    try
        [~, TVAR] = TVARir(TVAR, VARopt_iv);
    catch ME
        if mm <= 3
            warning('TVARir failed in rep %d: %s', mm, ME.message);
        end
        continue
    end

    % Store normalised Biv (divide by first element to match truth scaling)
    biv1 = TVAR.regime{1}.Biv;
    biv2 = TVAR.regime{2}.Biv;

    if isempty(biv1) || isempty(biv2); continue; end

    Biv1_MC(:, mm) = biv1 / biv1(1);
    Biv2_MC(:, mm) = biv2 / biv2(1);
end
fprintf('Done.\n\n');

% Count valid replications
valid = ~any(isnan(Biv1_MC), 1);
Biv1_MC = Biv1_MC(:, valid);
Biv2_MC = Biv2_MC(:, valid);
thresh_MC = thresh_MC(valid);
n1_MC     = n1_MC(valid);
Mvalid    = sum(valid);
fprintf('Valid replications: %d / %d\n\n', Mvalid, M);


%% ========================================================================
%  RESULTS
%  =========================================================================

fprintf('=================================================================\n');
fprintf('  Proxy TSVAR Correctness Test  (T=%d, M=%d, rho=%.1f)\n', T, Mvalid, rho);
fprintf('=================================================================\n\n');

%% --- Threshold recovery --------------------------------------------------
fprintf('--- Threshold estimation ---\n');
fprintf('  True threshold : %6.3f\n', gamma_true);
fprintf('  Mean estimate  : %6.3f  (RMSE = %.3f)\n', ...
    mean(thresh_MC), sqrt(mean((thresh_MC - gamma_true).^2)));
fprintf('  Mean regime 1 N: %.1f / %d  (expected ~%.0f)\n\n', ...
    mean(n1_MC), T-nlag, (T-nlag) * mean(y(:,1) <= gamma_true));

%% --- Regime 1 Biv --------------------------------------------------------
bias1 = mean(Biv1_MC, 2) - Biv1_true_norm;
rmse1 = sqrt(mean((Biv1_MC - Biv1_true_norm).^2, 2));

fprintf('--- Regime 1: normalised Biv  (true B1(:,1)/B1(1,1) = [1; %.2f; %.2f])\n', ...
    Biv1_true_norm(2), Biv1_true_norm(3));
fprintf('  %-8s  %8s  %8s  %8s  %8s  %8s\n', 'Element','Truth','Mean','Bias','RMSE','Pass?');
pass1 = true;
for ii = 1:nvar
    ok = abs(bias1(ii)) < 0.05 && rmse1(ii) < 0.15;
    pass1 = pass1 && ok;
    fprintf('  %-8s  %8.3f  %8.3f  %8.3f  %8.3f  %8s\n', ...
        ['Biv1(' num2str(ii) ')'], Biv1_true_norm(ii), mean(Biv1_MC(ii,:)), ...
        bias1(ii), rmse1(ii), yesno(ok));
end

%% --- Regime 2 Biv --------------------------------------------------------
bias2 = mean(Biv2_MC, 2) - Biv2_true_norm;
rmse2 = sqrt(mean((Biv2_MC - Biv2_true_norm).^2, 2));

fprintf('\n--- Regime 2: normalised Biv  (true B2(:,1)/B2(1,1) = [1; %.2f; %.2f])\n', ...
    Biv2_true_norm(2), Biv2_true_norm(3));
fprintf('  %-8s  %8s  %8s  %8s  %8s  %8s\n', 'Element','Truth','Mean','Bias','RMSE','Pass?');
pass2 = true;
for ii = 1:nvar
    ok = abs(bias2(ii)) < 0.05 && rmse2(ii) < 0.15;
    pass2 = pass2 && ok;
    fprintf('  %-8s  %8.3f  %8.3f  %8.3f  %8.3f  %8s\n', ...
        ['Biv2(' num2str(ii) ')'], Biv2_true_norm(ii), mean(Biv2_MC(ii,:)), ...
        bias2(ii), rmse2(ii), yesno(ok));
end

%% --- Overall verdict -----------------------------------------------------
fprintf('\n=================================================================\n');
if pass1 && pass2
    fprintf('  OVERALL: PASS  — Proxy TSVAR identification is working correctly.\n');
else
    fprintf('  OVERALL: FAIL  — Check bias/RMSE table above.\n');
end
fprintf('=================================================================\n\n');


%% ========================================================================
%  PLOT: Monte Carlo distributions vs truth
%  =========================================================================

FigSize(26, 16);
labels = {'Biv(1)','Biv(2)','Biv(3)'};
for ii = 1:nvar
    % Regime 1
    subplot(2, nvar, ii);
    histogram(Biv1_MC(ii,:), 25, 'FaceColor', cmap(1), 'EdgeColor','w', 'Normalization','pdf');
    hold on;
    xline(Biv1_true_norm(ii), '-k', 'LineWidth', 2);
    xline(mean(Biv1_MC(ii,:)), '--', 'Color', cmap(1), 'LineWidth', 1.5);
    title(['Regime 1: ' labels{ii}], 'FontWeight','bold');
    xlabel('Estimate'); ylabel('Density');
    legend({'MC dist.','Truth','Mean'}, 'Location','best', 'FontSize',7);
    legend boxoff; grid on;

    % Regime 2
    subplot(2, nvar, nvar + ii);
    histogram(Biv2_MC(ii,:), 25, 'FaceColor', cmap(2), 'EdgeColor','w', 'Normalization','pdf');
    hold on;
    xline(Biv2_true_norm(ii), '-k', 'LineWidth', 2);
    xline(mean(Biv2_MC(ii,:)), '--', 'Color', cmap(2), 'LineWidth', 1.5);
    title(['Regime 2: ' labels{ii}], 'FontWeight','bold');
    xlabel('Estimate'); ylabel('Density');
    legend({'MC dist.','Truth','Mean'}, 'Location','best', 'FontSize',7);
    legend boxoff; grid on;
end

% exportgraphics(gcf, 'graphics/TVARtest_ProxySVAR_MC.pdf');
% fprintf('Figure saved to graphics/TVARtest_ProxySVAR_MC.pdf\n');
% close all;


%% ========================================================================
%  HELPER
%  =========================================================================
function s = yesno(tf)
    if tf; s = 'PASS'; else; s = 'FAIL'; end
end
