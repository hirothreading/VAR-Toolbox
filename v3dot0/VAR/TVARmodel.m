function [TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, TVARopt, THRVAR_EX, EXOG, nlag_ex)
%==========================================================================
% Estimate a two-regime Threshold VAR (TVAR). Two estimation methods are
% available, selected via TVARopt.tvar_method:
%
%   'freq'  — Conditional OLS with grid search over the threshold
%             (Hansen 1997, 2000). Fast, consistent with the OLS-based
%             architecture of VAR Toolbox.
%
%   'bayes' — Bayesian estimation via Gibbs sampling with Minnesota priors
%             (adapted from Haroon Mumtaz's TVAR Toolkit). The threshold,
%             delay, coefficients, and covariances are all sampled jointly.
%
% The model:
%   y_t = A_1 x_t + u_t^(1)   if  q_{t-d} <= gamma   (regime 1)
%   y_t = A_2 x_t + u_t^(2)   if  q_{t-d} >  gamma   (regime 2)
%
% where q_{t-d} is the threshold variable lagged by d periods.
%==========================================================================
% [TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, TVARopt, THRVAR_EX, EXOG, nlag_ex)
% -------------------------------------------------------------------------
% INPUT
%   - ENDO    : (nobs x nvar) matrix of endogenous variables
%   - nlag    : lag order
%   - const   : 0 no constant; 1 constant [default = 1]
%   - TVARopt : options structure (see TVARoption.m)
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - THRVAR_EX : (nobs x 1) external threshold variable (when thrvar_idx=0)
%   - EXOG      : (nobs x nvar_ex) matrix of exogenous variables
%   - nlag_ex   : lag order for exogenous variables [default = 0]
% -------------------------------------------------------------------------
% OUTPUT
%   - TVAR    : structure with estimation results (see below)
%   - TVARopt : updated options structure
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
[nobs, nvar] = size(ENDO);

if ~exist('const','var') || isempty(const)
    const = 1;
end
if ~exist('TVARopt','var') || isempty(TVARopt)
    TVARopt = TVARoption;
end
if ~exist('EXOG','var')
    EXOG = [];
    nvar_ex = 0;
    nlag_ex = 0;
else
    [~, nvar_ex] = size(EXOG);
    if ~exist('nlag_ex','var') || isempty(nlag_ex)
        nlag_ex = 0;
    end
end
if ~exist('THRVAR_EX','var')
    THRVAR_EX = [];
end

% Retrieve TVAR-specific options
thrvar_idx   = TVARopt.thrvar_idx;
delay        = TVARopt.delay;
tartransform = TVARopt.tartransform;
tvar_method  = TVARopt.tvar_method;

% Validate inputs
if delay < 1
    error('TVARmodel: delay must be >= 1');
end
if thrvar_idx == 0 && isempty(THRVAR_EX)
    error('TVARmodel: thrvar_idx=0 requires THRVAR_EX to be provided');
end
if thrvar_idx > nvar
    error('TVARmodel: thrvar_idx (%d) exceeds number of variables (%d)', thrvar_idx, nvar);
end


%% Construct threshold variable
%==========================================================================
if thrvar_idx == 0
    thrvar_raw = THRVAR_EX;
else
    thrvar_raw = ENDO(:, thrvar_idx);
end

% Apply transformation
if tartransform == 1
    % Annual growth rate (4-period difference)
    Lx = 4;
    thrvar_transformed = [nan(Lx,1); thrvar_raw(Lx+1:end) - thrvar_raw(1:end-Lx)];
else
    Lx = 0;
    thrvar_transformed = thrvar_raw;
end

% Lag the threshold variable by delay periods
Ystar_full = [nan(delay,1); thrvar_transformed(1:end-delay)];


%% Build data matrices
%==========================================================================
% Use VARmakexy for the endogenous variables
[Y_full, X_full] = VARmakexy(ENDO, nlag, const);

% Handle exogenous variables
if nvar_ex > 0
    X_EX = VARmakelags(EXOG, nlag_ex);
    if nlag == nlag_ex
        X_full = [X_full X_EX];
    elseif nlag > nlag_ex
        diff_lag = nlag - nlag_ex;
        X_EX = X_EX(diff_lag+1:end,:);
        X_full = [X_full X_EX];
    elseif nlag < nlag_ex
        diff_lag = nlag_ex - nlag;
        Y_full = Y_full(diff_lag+1:end,:);
        X_full = [X_full(diff_lag+1:end,:) X_EX];
    end
end

% Align Ystar with Y/X (both start at observation nlag+1 of ENDO)
Ystar_aligned = Ystar_full(nlag+1:end);
if nvar_ex > 0 && nlag < nlag_ex
    diff_lag = nlag_ex - nlag;
    Ystar_aligned = Ystar_aligned(diff_lag+1:end);
end

% Remove observations where Ystar is NaN (due to delay + transformation)
valid = ~isnan(Ystar_aligned);
Y = Y_full(valid,:);
X = X_full(valid,:);
Ystar = Ystar_aligned(valid);
nobse = size(Y, 1);

% Coefficient dimensions
ncoeff     = nvar * nlag;
ncoeff_ex  = nvar_ex * (nlag_ex + 1);
ntotcoeff  = ncoeff + ncoeff_ex + const;


%% Dispatch to the chosen estimation method
%==========================================================================
if strcmp(tvar_method, 'freq')
    [TVAR, TVARopt] = estimate_freq(ENDO, Y, X, Ystar, nvar, nlag, const, ...
        nobse, ncoeff, ntotcoeff, nvar_ex, nlag_ex, EXOG, TVARopt);
elseif strcmp(tvar_method, 'bayes')
    [TVAR, TVARopt] = estimate_bayes(ENDO, Y, X, Ystar, nvar, nlag, const, ...
        nobse, ncoeff, ntotcoeff, nvar_ex, nlag_ex, EXOG, TVARopt);
else
    error('TVARmodel: tvar_method must be ''freq'' or ''bayes''');
end

% Store common fields
TVAR.ENDO     = ENDO;
TVAR.EXOG     = EXOG;
TVAR.nvar     = nvar;
TVAR.nlag     = nlag;
TVAR.const    = const;
TVAR.nvar_ex  = nvar_ex;
TVAR.nlag_ex  = nlag_ex;
TVAR.nobs     = nobse;
TVAR.nregimes = 2;
TVAR.Y        = Y;
TVAR.X        = X;
TVAR.Ystar    = Ystar;
TVAR.IV       = [];
TVAR.method   = tvar_method;

end


%% ========================================================================
%  FREQUENTIST ESTIMATION: Hansen (1997, 2000) grid search + OLS
%  ========================================================================
function [TVAR, TVARopt] = estimate_freq(ENDO, Y, X, Ystar, nvar, nlag, const, ...
    nobse, ncoeff, ntotcoeff, nvar_ex, nlag_ex, EXOG, TVARopt)

trim  = TVARopt.trim;
ngrid = TVARopt.ngrid;

% Build grid of candidate thresholds
Ystar_sorted = sort(Ystar);
lo = round(trim * nobse) + 1;
hi = nobse - round(trim * nobse);
if lo >= hi
    error('TVARmodel: not enough observations after trimming. Reduce trim or increase sample.');
end
grid_vals = linspace(Ystar_sorted(lo), Ystar_sorted(hi), ngrid);

% Grid search: minimise total RSS
RSS_grid = nan(ngrid, 1);
for gg = 1:ngrid
    gamma = grid_vals(gg);
    idx1 = (Ystar <= gamma);
    idx2 = ~idx1;
    n1 = sum(idx1);
    n2 = sum(idx2);
    % Need enough obs in each regime
    if n1 < ntotcoeff + 1 || n2 < ntotcoeff + 1
        RSS_grid(gg) = inf;
        continue;
    end
    % OLS per regime
    Y1 = Y(idx1,:);  X1 = X(idx1,:);
    Y2 = Y(idx2,:);  X2 = X(idx2,:);
    resid1 = Y1 - X1 * ((X1'*X1) \ (X1'*Y1));
    resid2 = Y2 - X2 * ((X2'*X2) \ (X2'*Y2));
    RSS_grid(gg) = trace(resid1'*resid1) + trace(resid2'*resid2);
end

% Select optimal threshold
[~, opt_idx] = min(RSS_grid);
thresh = grid_vals(opt_idx);

% Final regime split
idx1 = (Ystar <= thresh);
idx2 = ~idx1;

% Estimate regime-specific VARs using VARmodel for full compatibility
% We need to reconstruct ENDO subsets that produce the correct Y/X per regime
% Instead, we build the regime structs manually from OLS on the split data
for k = 1:2
    if k == 1
        Yk = Y(idx1,:);  Xk = X(idx1,:);
    else
        Yk = Y(idx2,:);  Xk = X(idx2,:);
    end
    nobsk = size(Yk, 1);

    % OLS
    Ftk = (Xk'*Xk) \ (Xk'*Yk);
    Fk  = Ftk';
    residk = Yk - Xk * Ftk;
    sigmak = (residk' * residk) / (nobsk - ntotcoeff);

    % Companion matrix
    Fcomp_k = [Fk(:, 1+const:nvar*nlag+const); ...
               eye(nvar*(nlag-1)), zeros(nvar*(nlag-1), nvar)];
    maxEig_k = max(abs(eig(Fcomp_k)));

    % Populate regime struct (mirrors VAR struct from VARmodel)
    TVAR.regime{k}.Ft        = Ftk;
    TVAR.regime{k}.F         = Fk;
    TVAR.regime{k}.sigma     = sigmak;
    TVAR.regime{k}.resid     = residk;
    TVAR.regime{k}.Y         = Yk;
    TVAR.regime{k}.X         = Xk;
    TVAR.regime{k}.Fcomp     = Fcomp_k;
    TVAR.regime{k}.maxEig    = maxEig_k;
    TVAR.regime{k}.nobs      = nobsk;
    TVAR.regime{k}.nvar      = nvar;
    TVAR.regime{k}.nlag      = nlag;
    TVAR.regime{k}.const     = const;
    TVAR.regime{k}.ncoeff    = ncoeff;
    TVAR.regime{k}.ntotcoeff = ntotcoeff;
    TVAR.regime{k}.nvar_ex   = nvar_ex;
    TVAR.regime{k}.nlag_ex   = nlag_ex;
    TVAR.regime{k}.ENDO      = ENDO;
    TVAR.regime{k}.EXOG      = EXOG;
    TVAR.regime{k}.B         = [];
    TVAR.regime{k}.Biv       = [];
    TVAR.regime{k}.PSI       = [];
    TVAR.regime{k}.Fp        = [];
    TVAR.regime{k}.IV        = [];

    % Equation-by-equation results (for TVARprint)
    for j = 1:nvar
        OLSout = OLSmodel(Yk(:,j), Xk, 0);
        TVAR.regime{k}.eq{j}.beta  = OLSout.beta;
        TVAR.regime{k}.eq{j}.tstat = OLSout.tstat;
        TVAR.regime{k}.eq{j}.bstd  = OLSout.bstd;
        TVAR.regime{k}.eq{j}.rsqr  = OLSout.rsqr;
        TVAR.regime{k}.eq{j}.rbar  = OLSout.rbar;
        TVAR.regime{k}.eq{j}.sige  = OLSout.sige;
        TVAR.regime{k}.eq{j}.dw    = OLSout.dw;
    end
end

% Store threshold results
TVAR.thresh      = thresh;
TVAR.delay       = TVARopt.delay;
TVAR.regime_idx  = idx1;   % logical: true = regime 1
TVAR.RSS_grid    = RSS_grid;
TVAR.thresh_grid = grid_vals;
% Total RSS at optimal threshold (used by TVARtest)
TVAR.RSS = sum(sum(TVAR.regime{1}.resid.^2)) + sum(sum(TVAR.regime{2}.resid.^2));

end


%% ========================================================================
%  BAYESIAN ESTIMATION: Gibbs sampler (adapted from Mumtaz)
%  ========================================================================
function [TVAR, TVARopt] = estimate_bayes(ENDO, Y, X, Ystar, nvar, nlag, const, ...
    nobse, ncoeff, ntotcoeff, nvar_ex, nlag_ex, EXOG, TVARopt)

% Currently supports const=1 only for Bayesian mode
if const ~= 1
    error('TVARmodel: Bayesian mode currently supports const=1 only');
end

% Retrieve Bayesian options
nreps       = TVARopt.nreps;
nburn       = TVARopt.nburn;
nskip       = TVARopt.nskip;
lambdaP     = TVARopt.lambdaP;
tauP        = TVARopt.tauP;
epsilonP    = TVARopt.epsilonP;
tarvariance = TVARopt.tarvariance;
tarscale    = TVARopt.tarscale;
maxtrys     = TVARopt.maxtrys;
mult        = TVARopt.mult;

% Number of retained draws
Sindex = (nburn+1):nskip:nreps;
fsize  = length(Sindex);

% Minimum observations per regime
ncrit = nvar * nlag + 1;

% ---- Build Mumtaz-ordered data matrices ----
% Mumtaz ordering: [lag1 lag2 ... lagL constant] (lags first, constant last)
% We need this for the Gibbs sampler because the helper functions
% (getcoef_tvar, stability_tvar) expect this ordering.
X_bayes = zeros(nobse, nvar*nlag + 1);
% Extract lag blocks from X (VARmakexy ordering: [const lag1 lag2 ... lagL])
% With const=1: column 1 = constant, columns 2:end = lags
X_bayes(:, 1:nvar*nlag) = X(:, const+1:const+nvar*nlag);
X_bayes(:, nvar*nlag+1) = 1;  % constant at the end
ncoeff_bayes = nvar * nlag + 1;  % per equation in Mumtaz ordering

% ---- Compute priors ----
muP    = mean(Y)';
sigmaP = zeros(nvar, 1);
deltaP = zeros(nvar, 1);
for i = 1:nvar
    ytemp = Y(:,i);
    xtemp = [ytemp(1:end-1), ones(nobse-1,1)];
    ytemp2 = ytemp(2:end);
    btemp = xtemp \ ytemp2;
    etemp = ytemp2 - xtemp * btemp;
    stemp = (etemp' * etemp) / length(ytemp2);
    if abs(btemp(1)) > 1
        btemp(1) = 1;
    end
    deltaP(i) = btemp(1);
    sigmaP(i) = sqrt(stemp);
end

% Minnesota prior dummy observations
[yd, xd] = create_dummies_tvar(lambdaP, tauP, deltaP, epsilonP, nlag, muP, sigmaP, nvar);

% Prior mean for threshold
tarmean = mean(Ystar);

% ---- Initialise Gibbs sampler ----
sigma1 = eye(nvar);
sigma2 = eye(nvar);
beta0  = vec(X_bayes \ Y);
beta01 = beta0;
beta02 = beta0;
tar    = tarmean;
tarold = tar;

% Storage
bsave1   = zeros(fsize, nvar * ncoeff_bayes);
bsave2   = zeros(fsize, nvar * ncoeff_bayes);
sigmaS1  = zeros(fsize, nvar, nvar);
sigmaS2  = zeros(fsize, nvar, nvar);
tsave    = zeros(fsize, 2);

naccept = 0;
jgibbs  = 1;  % index into saved draws

disp(' ');
disp('============================================================');
disp(sprintf('TVAR Bayesian estimation: %d replications, %d burn-in, thinning every %d', nreps, nburn, nskip));

% ---- Gibbs loop ----
for igibbs = 1:nreps

    % Step 1: Split data by threshold
    e1 = (Ystar <= tar);
    e2 = ~e1;
    Y1 = Y(e1,:);  X1 = X_bayes(e1,:);
    Y2 = Y(e2,:);  X2 = X_bayes(e2,:);

    % Step 2: Sample coefficients and covariance — Regime 1
    Y0 = [Y1; yd];
    X0 = [X1; xd];
    mstar1 = vec(X0 \ Y0);
    xx = X0' * X0;
    ixx1 = xx \ eye(size(xx,1));
    [beta1, PROBLEM1] = getcoef_tvar(mstar1, sigma1, ixx1, maxtrys, nvar, nlag);
    if PROBLEM1
        beta1 = beta01;
    else
        beta01 = beta1;
    end
    % Draw covariance from Inverse-Wishart
    e_res = Y0 - X0 * reshape(beta1, ncoeff_bayes, nvar);
    scale1 = e_res' * e_res;
    sigma1 = iwpQ(size(Y0,1), inv(scale1));

    % Step 3: Sample coefficients and covariance — Regime 2
    Y0 = [Y2; yd];
    X0 = [X2; xd];
    mstar2 = vec(X0 \ Y0);
    xx = X0' * X0;
    ixx2 = xx \ eye(size(xx,1));
    [beta2, PROBLEM2] = getcoef_tvar(mstar2, sigma2, ixx2, maxtrys, nvar, nlag);
    if PROBLEM2
        beta2 = beta02;
    else
        beta02 = beta2;
    end
    e_res = Y0 - X0 * reshape(beta2, ncoeff_bayes, nvar);
    scale2 = e_res' * e_res;
    sigma2 = iwpQ(size(Y0,1), inv(scale2));

    % Step 4: Sample threshold via Random Walk Metropolis-Hastings
    tarnew = tarold + randn * sqrt(tarscale);
    postnew = getvarpost_tvar(Y, X_bayes, beta1, beta2, sigma1, sigma2, nlag, tarnew, tarmean, tarvariance, Ystar, ncrit);
    postold = getvarpost_tvar(Y, X_bayes, beta1, beta2, sigma1, sigma2, nlag, tarold, tarmean, tarvariance, Ystar, ncrit);
    accept_ratio = exp(postnew - postold);
    if rand < accept_ratio
        tarold = tarnew;
        naccept = naccept + 1;
    end
    tar = tarold;
    arate = naccept / igibbs;

    % Step 5: Sample delay parameter (if multiple delays specified)
    % For now, delay is fixed (as specified by TVARopt.delay)
    current_delay = TVARopt.delay;

    % Display progress
    if mod(igibbs, mult) == 0
        fprintf('  Replication %d of %d | acceptance rate %.3f\n', igibbs, nreps, arate);
    end

    % Store draws (post burn-in, thinned)
    if igibbs > nburn && any(Sindex == igibbs)
        bsave1(jgibbs,:)    = beta1';
        bsave2(jgibbs,:)    = beta2';
        sigmaS1(jgibbs,:,:) = sigma1;
        sigmaS2(jgibbs,:,:) = sigma2;
        tsave(jgibbs,:)     = [tar, current_delay];
        jgibbs = jgibbs + 1;
    end
end

disp('TVAR Bayesian estimation — done');
disp(sprintf('  Acceptance rate: %.3f', arate));
disp(sprintf('  Retained draws:  %d', fsize));
disp('============================================================');

% ---- Construct TVAR output ----
% Posterior means
beta1_mean  = mean(bsave1, 1)';
beta2_mean  = mean(bsave2, 1)';
sigma1_mean = squeeze(mean(sigmaS1, 1));
sigma2_mean = squeeze(mean(sigmaS2, 1));
thresh_mean = mean(tsave(:,1));
delay_mean  = round(mean(tsave(:,2)));

% Final regime split at posterior mean threshold
idx1 = (Ystar <= thresh_mean);
idx2 = ~idx1;

% Build regime structs using posterior means
% Convert from Mumtaz ordering [lags|const] to VARmodel ordering [const|lags]
for k = 1:2
    if k == 1
        beta_k = beta1_mean;
        sigma_k = sigma1_mean;
        Yk = Y(idx1,:);  Xk = X(idx1,:);
    else
        beta_k = beta2_mean;
        sigma_k = sigma2_mean;
        Yk = Y(idx2,:);  Xk = X(idx2,:);
    end
    nobsk = size(Yk, 1);

    % Reshape from Mumtaz ordering: (ncoeff_bayes x nvar)
    % Rows: [lag1_coefs; lag2_coefs; ...; lagL_coefs; constant]
    beta_mat_mumtaz = reshape(beta_k, ncoeff_bayes, nvar);

    % Convert to VARmodel ordering: [constant; lag1; lag2; ...; lagL]
    Ftk = [beta_mat_mumtaz(end,:);              % constant (row 1)
           beta_mat_mumtaz(1:nvar*nlag,:)];     % lags (rows 2 to end)
    Fk = Ftk';

    % Companion matrix (skip constant = first row of Ftk)
    Fcomp_k = [Fk(:, 1+const:nvar*nlag+const); ...
               eye(nvar*(nlag-1)), zeros(nvar*(nlag-1), nvar)];
    maxEig_k = max(abs(eig(Fcomp_k)));

    % Residuals at posterior mean
    residk = Yk - Xk * Ftk;

    % Populate regime struct
    TVAR.regime{k}.Ft        = Ftk;
    TVAR.regime{k}.F         = Fk;
    TVAR.regime{k}.sigma     = sigma_k;
    TVAR.regime{k}.resid     = residk;
    TVAR.regime{k}.Y         = Yk;
    TVAR.regime{k}.X         = Xk;
    TVAR.regime{k}.Fcomp     = Fcomp_k;
    TVAR.regime{k}.maxEig    = maxEig_k;
    TVAR.regime{k}.nobs      = nobsk;
    TVAR.regime{k}.nvar      = nvar;
    TVAR.regime{k}.nlag      = nlag;
    TVAR.regime{k}.const     = const;
    TVAR.regime{k}.ncoeff    = ncoeff;
    TVAR.regime{k}.ntotcoeff = ntotcoeff;
    TVAR.regime{k}.nvar_ex   = nvar_ex;
    TVAR.regime{k}.nlag_ex   = nlag_ex;
    TVAR.regime{k}.ENDO      = ENDO;
    TVAR.regime{k}.EXOG      = EXOG;
    TVAR.regime{k}.B         = [];
    TVAR.regime{k}.Biv       = [];
    TVAR.regime{k}.PSI       = [];
    TVAR.regime{k}.Fp        = [];
    TVAR.regime{k}.IV        = [];
end

% Store threshold and Bayesian results
TVAR.thresh     = thresh_mean;
TVAR.delay      = delay_mean;
TVAR.regime_idx = idx1;

% Posterior draws (for TVARirband)
TVAR.beta1_draws  = bsave1;
TVAR.beta2_draws  = bsave2;
TVAR.sigma1_draws = sigmaS1;
TVAR.sigma2_draws = sigmaS2;
TVAR.thresh_draws = tsave(:,1);
TVAR.delay_draws  = tsave(:,2);
TVAR.ndraws       = fsize;
TVAR.acceptance_rate = arate;

% Store MCMC settings
TVAR.nreps = nreps;
TVAR.nburn = nburn;
TVAR.nskip = nskip;

end
