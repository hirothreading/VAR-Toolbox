function [TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, THRVAR_EX, EXOG, nlag_ex)
%==========================================================================
% Estimate a two-regime Threshold VAR (TVAR) by conditional OLS with grid
% search over the threshold value (Hansen 1997, 2000).
%
% The model switches between two regimes based on whether q_{t-d} is below
% or above the estimated threshold gamma:
%
%   y_t = A_1 x_t  +  u_t^(1)   if  q_{t-d} <= gamma   (regime 1)
%   y_t = A_2 x_t  +  u_t^(2)   if  q_{t-d} >  gamma   (regime 2)
%
% where q_{t-d} is the threshold variable lagged by d periods (the delay).
%==========================================================================
% [TVAR, TVARopt] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, THRVAR_EX, EXOG, nlag_ex)
% -------------------------------------------------------------------------
% INPUT
%   - ENDO       : (nobs x nvar) matrix of endogenous variables
%   - nlag       : lag order of the VAR
%   - const      : 0 no constant; 1 constant; 2 constant+trend;
%                  3 constant+trend+trend^2  [default = 1]
%   - thrvar_idx : column index in ENDO used as threshold variable.
%                  Set to 0 if providing an external threshold via THRVAR_EX.
%   - delay      : delay d for the threshold variable; q_{t-d} is used.
%                  Must satisfy 1 <= delay <= nlag.  [default = 1]
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - THRVAR_EX  : (nobs x 1) external threshold variable; required when
%                  thrvar_idx = 0, ignored otherwise.
%   - EXOG       : (nobs x nvar_ex) matrix of exogenous variables
%   - nlag_ex    : lag order for exogenous variables  [default = 0]
% -------------------------------------------------------------------------
% OUTPUT
%   - TVAR       : structure with TVAR estimation results. The sub-struct
%                  TVAR.regime{k} (k=1,2) is shaped identically to the VAR
%                  struct returned by VARmodel, making it compatible with
%                  VARir, VARvd, and related functions.
%   - TVARopt    : structure with TVAR options (see TVARoption)
% -------------------------------------------------------------------------
% EXAMPLE
%   - See TVARToolbox_Primer.m in "../Primer/"
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================


%% Check inputs
%--------------------------------------------------------------------------
[nobs, nvar] = size(ENDO);

% Create TVARopt
TVARopt = TVARoption;

% const
if ~exist('const','var') || isempty(const)
    const = 1;
end

% thrvar_idx
if ~exist('thrvar_idx','var') || isempty(thrvar_idx)
    thrvar_idx = TVARopt.thrvar_idx;
end
TVARopt.thrvar_idx = thrvar_idx;

% delay
if ~exist('delay','var') || isempty(delay)
    delay = TVARopt.delay;
end
if delay < 1
    error('TVARmodel: delay must be >= 1');
end
TVARopt.delay = delay;

% External threshold variable
if thrvar_idx == 0
    if ~exist('THRVAR_EX','var') || isempty(THRVAR_EX)
        error('TVARmodel: thrvar_idx=0 requires THRVAR_EX to be provided');
    end
    if size(THRVAR_EX,1) ~= nobs
        error('TVARmodel: THRVAR_EX must have the same number of rows as ENDO');
    end
    use_external_thr = true;
else
    if thrvar_idx < 1 || thrvar_idx > nvar
        error('TVARmodel: thrvar_idx must be between 1 and nvar (or 0 for external)');
    end
    use_external_thr = false;
end

% Exogenous variables
if exist('EXOG','var') && ~isempty(EXOG)
    [nobs2, nvar_ex] = size(EXOG);
    if nobs2 ~= nobs
        error('TVARmodel: nobs in EXOG-matrix not the same as ENDO');
    end
    if ~exist('nlag_ex','var') || isempty(nlag_ex)
        nlag_ex = 0;
    end
    has_exog = true;
else
    nvar_ex  = 0;
    nlag_ex  = 0;
    EXOG     = [];
    has_exog = false;
end


%% Build Y and X matrices (identical to VARmodel)
%--------------------------------------------------------------------------
nobse        = nobs - max(nlag, nlag_ex);
ncoeff       = nvar * nlag;
ncoeff_ex    = nvar_ex * (nlag_ex + 1);
ntotcoeff    = ncoeff + ncoeff_ex + const;

[Y, X] = VARmakexy(ENDO, nlag, const);

if has_exog
    X_EX = VARmakelags(EXOG, nlag_ex);
    if nlag == nlag_ex
        X = [X X_EX];
    elseif nlag > nlag_ex
        diff_lags = nlag - nlag_ex;
        X_EX = X_EX(diff_lags+1:end, :);
        X = [X X_EX];
    else
        diff_lags = nlag_ex - nlag;
        Y = Y(diff_lags+1:end, :);
        X = [X(diff_lags+1:end,:) X_EX];
    end
end

% Recompute nobse after possible lag adjustment
nobse = size(Y, 1);


%% Build threshold variable
%--------------------------------------------------------------------------
% The threshold variable for the t-th row of Y (which corresponds to
% calendar time nlag+t in ENDO) is q_{nlag+t-delay} = ENDO(nlag+t-delay).
% In vector form: rows nlag+1-delay through nobs-delay of the threshold col.
thr_start = nlag + 1 - delay;
thr_end   = nobs - delay;

if thr_start < 1
    error('TVARmodel: delay (%d) exceeds nlag (%d). Reduce delay or increase nlag.', delay, nlag);
end

if use_external_thr
    q = THRVAR_EX(thr_start:thr_end);
else
    q = ENDO(thr_start:thr_end, thrvar_idx);
end
q = q(:);  % ensure column vector

if length(q) ~= nobse
    error('TVARmodel: threshold variable length (%d) does not match nobse (%d)', length(q), nobse);
end


%% Grid search for optimal threshold
%--------------------------------------------------------------------------
trim  = TVARopt.trim;
ngrid = TVARopt.ngrid;

q_lo = quantile(q, trim);
q_hi = quantile(q, 1 - trim);

if q_lo >= q_hi
    error('TVARmodel: threshold grid is empty (trim too aggressive or too little variation in q)');
end

grid_gamma = linspace(q_lo, q_hi, ngrid);
grid_RSS   = nan(ngrid, 1);

min_obs_per_regime = ntotcoeff + 1;  % minimum observations for OLS

for jj = 1:ngrid
    gamma_j = grid_gamma(jj);
    idx1 = find(q <= gamma_j);
    idx2 = find(q >  gamma_j);

    % Skip if either regime is too small to identify
    if length(idx1) < min_obs_per_regime || length(idx2) < min_obs_per_regime
        continue
    end

    X1 = X(idx1, :);  Y1 = Y(idx1, :);
    X2 = X(idx2, :);  Y2 = Y(idx2, :);

    Ft1 = (X1'*X1) \ (X1'*Y1);
    Ft2 = (X2'*X2) \ (X2'*Y2);

    resid1 = Y1 - X1*Ft1;
    resid2 = Y2 - X2*Ft2;

    grid_RSS(jj) = sum(sum(resid1.^2)) + sum(sum(resid2.^2));
end

if all(isnan(grid_RSS))
    error('TVARmodel: no valid threshold found — all grid points have insufficient regime observations. Try reducing trim or increasing the sample.');
end

% Optimal threshold: minimise total RSS
[RSS_opt, idx_opt] = min(grid_RSS);
thresh = grid_gamma(idx_opt);

% Final regime assignments
regime_idx      = ones(nobse, 1);
regime_idx(q > thresh) = 2;
obs_idx1 = find(regime_idx == 1);
obs_idx2 = find(regime_idx == 2);


%% Final OLS estimation by regime
%--------------------------------------------------------------------------
for k = 1:2
    if k == 1
        obs_idx = obs_idx1;
    else
        obs_idx = obs_idx2;
    end

    Xk = X(obs_idx, :);
    Yk = Y(obs_idx, :);
    nobs_k = length(obs_idx);

    % System-level OLS
    Ft_k = (Xk'*Xk) \ (Xk'*Yk);
    F_k  = Ft_k';
    resid_k = Yk - Xk*Ft_k;
    sigma_k = (1/(nobs_k - ntotcoeff)) * (resid_k' * resid_k);

    % Companion matrix
    Fcomp_k = [F_k(:, 1+const : nvar*nlag+const); ...
               eye(nvar*(nlag-1)) zeros(nvar*(nlag-1), nvar)];
    maxEig_k = max(abs(eig(Fcomp_k)));

    % Equation-by-equation OLS (for TVARprint and eq1..eqN fields)
    reg = struct();
    for j = 1:nvar
        Yvec   = Yk(:, j);
        OLSout = OLSmodel(Yvec, Xk, 0);
        eqname = ['eq' num2str(j)];
        reg.(eqname).beta  = OLSout.beta;
        reg.(eqname).tstat = OLSout.tstat;
        reg.(eqname).bstd  = OLSout.bstd;
        reg.(eqname).tprob = OLSout.tprob;
        reg.(eqname).resid = OLSout.resid;
        reg.(eqname).yhat  = OLSout.yhat;
        reg.(eqname).y     = Yvec;
        reg.(eqname).rsqr  = OLSout.rsqr;
        reg.(eqname).rbar  = OLSout.rbar;
        reg.(eqname).sige  = OLSout.sige;
        reg.(eqname).dw    = OLSout.dw;
    end

    % Assemble regime sub-struct (VARir-compatible)
    reg.Ft         = Ft_k;
    reg.F          = F_k;
    reg.sigma      = sigma_k;
    reg.resid      = resid_k;
    reg.X          = Xk;
    reg.Y          = Yk;
    reg.obs_idx    = obs_idx;   % row indices into full Y/X (used by TVARir for iv)
    reg.Fcomp      = Fcomp_k;
    reg.maxEig     = maxEig_k;
    reg.nobs       = nobs_k;
    reg.nvar       = nvar;
    reg.nlag       = nlag;
    reg.const      = const;
    reg.ntotcoeff  = ntotcoeff;
    reg.ncoeff     = ncoeff;
    reg.ncoeff_ex  = ncoeff_ex;
    reg.nvar_ex    = nvar_ex;
    reg.nlag_ex    = nlag_ex;
    reg.ENDO       = ENDO;      % full data (needed by TVARirband)
    reg.EXOG       = EXOG;
    % Placeholders populated later by TVARir
    reg.B          = [];
    reg.Biv        = [];
    reg.PSI        = [];
    reg.Fp         = [];
    reg.IV         = [];

    TVAR.regime{k} = reg;
end


%% Assemble TVAR struct
%--------------------------------------------------------------------------
TVAR.nthresh     = 1;
TVAR.nregimes    = 2;
TVAR.thresh      = thresh;
TVAR.delay       = delay;
TVAR.thrvar_idx  = thrvar_idx;
TVAR.thrvar      = q;           % (nobse x 1) threshold variable
TVAR.regime_idx  = regime_idx;  % (nobse x 1) regime assignment (1 or 2)
TVAR.nobs        = nobse;
TVAR.nvar        = nvar;
TVAR.nlag        = nlag;
TVAR.const       = const;
TVAR.ntotcoeff   = ntotcoeff;
TVAR.ENDO        = ENDO;
TVAR.EXOG        = EXOG;
TVAR.X           = X;
TVAR.Y           = Y;
TVAR.IV          = [];          % populated by user before TVARir (for iv identification)
TVAR.RSS         = RSS_opt;
TVAR.grid_gamma  = grid_gamma;
TVAR.grid_RSS    = grid_RSS;
if has_exog
    TVAR.X_EX    = X_EX;
end
