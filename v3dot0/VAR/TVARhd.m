function [HD, TVAR] = TVARhd(TVAR, VARopt)
%==========================================================================
% Compute the historical decomposition for a Threshold VAR estimated with
% TVARmodel. Structural impact matrices and companion matrices switch with
% the regime at each observation.
%==========================================================================
% [HD, TVAR] = TVARhd(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel
%   - VARopt : VAR options (ident, vnames, etc.)
% -------------------------------------------------------------------------
% OUTPUT
%   - HD     : structure with fields:
%       .shock  (nobs+nlag x nvar x nvar)   contribution of each shock
%       .init   (nobs+nlag x nvar)          initial-condition contribution
%       .const  (nobs+nlag x nvar)          constant contribution
%       .trend  (nobs+nlag x nvar)          linear trend contribution
%       .trend2 (nobs+nlag x nvar)          quadratic trend contribution
%       .endo   (nobs+nlag x nvar)          sum of all components (= data)
%   - TVAR   : updated with B matrices per regime
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('TVAR','var')
    error('You need to provide TVAR structure, result of TVARmodel');
end
if strcmp(VARopt.ident,'iv')
    disp('---------------------------------------------')
    disp('Historical decomposition not available with')
    disp('external instruments identification (iv)');
    disp('---------------------------------------------')
    error('ERROR. See details above');
end


%% Retrieve variables
%==========================================================================
nvar     = TVAR.nvar;
nlag     = TVAR.nlag;
const    = TVAR.const;
nvar_ex  = TVAR.nvar_ex;
nvarXeq  = nvar * nlag;
Y        = TVAR.Y;                              % (nobse x nvar)
X        = TVAR.X(:, 1+const:nvarXeq+const);    % lagged endogenous only
nobs     = size(Y, 1);
idx1     = TVAR.regime_idx;                      % logical (nobse x 1)
nregimes = TVAR.nregimes;


%% Identification: Recover B matrix per regime
%==========================================================================
B     = cell(nregimes, 1);
Fcomp = cell(nregimes, 1);
F     = cell(nregimes, 1);

for k = 1:nregimes
    sigma_k = TVAR.regime{k}.sigma;
    Fcomp{k} = TVAR.regime{k}.Fcomp;
    F{k}     = TVAR.regime{k}.Ft';    % make comparable to notes (nvar x ntotcoeff)

    if strcmp(VARopt.ident, 'short')
        [out, chol_flag] = chol(sigma_k);
        if chol_flag ~= 0
            error('VCV of regime %d is not positive definite', k);
        end
        B{k} = out';

    elseif strcmp(VARopt.ident, 'long')
        Finf_big = inv(eye(length(Fcomp{k})) - Fcomp{k});
        Finf = Finf_big(1:nvar, 1:nvar);
        D = chol(Finf * sigma_k * Finf')';
        B{k} = Finf \ D;

    elseif strcmp(VARopt.ident, 'sign')
        if isempty(TVAR.regime{k}.B)
            error('Regime %d: B matrix not set. Use SR.m first.', k);
        end
        B{k} = TVAR.regime{k}.B;

    else
        error('Identification incorrectly specified. Choose: short, long, sign.');
    end
end


%% Compute full residual series using regime-switching coefficients
%==========================================================================
resid_full = zeros(nobs, nvar);
for t = 1:nobs
    if idx1(t)
        resid_full(t,:) = Y(t,:) - TVAR.X(t,:) * TVAR.regime{1}.Ft;
    else
        resid_full(t,:) = Y(t,:) - TVAR.X(t,:) * TVAR.regime{2}.Ft;
    end
end


%% Regime index function: returns regime k for observation t
%==========================================================================
% k(t) = 1 if idx1(t)==true, 2 otherwise
kfun = @(t) 2 - idx1(t);  % idx1=true -> k=1; idx1=false -> k=2


%% Compute historical decompositions with regime switching
%==========================================================================
Icomp = [eye(nvar) zeros(nvar, (nlag-1)*nvar)];

% --- Contribution of each shock ---
HDshock_big = zeros(nlag*nvar, nobs+1, nvar);
HDshock     = zeros(nvar, nobs+1, nvar);
for t = 1:nobs
    k = kfun(t);
    eps_t = B{k} \ resid_full(t,:)';       % structural shocks at t
    for j = 1:nvar
        eps_big = zeros(nlag*nvar, 1);
        eps_big(j) = eps_t(j);
        B_big = zeros(nvarXeq, nvar);
        B_big(1:nvar,:) = B{k};
        HDshock_big(:,t+1,j) = B_big * eps_big(1:nvar) + Fcomp{k} * HDshock_big(:,t,j);
        HDshock(:,t+1,j) = Icomp * HDshock_big(:,t+1,j);
    end
end

% --- Initial value ---
HDinit_big = zeros(nlag*nvar, nobs+1);
HDinit     = zeros(nvar, nobs+1);
HDinit_big(:,1) = X(1,:)';
HDinit(:,1) = Icomp * HDinit_big(:,1);
for t = 1:nobs
    k = kfun(t);
    HDinit_big(:,t+1) = Fcomp{k} * HDinit_big(:,t);
    HDinit(:,t+1) = Icomp * HDinit_big(:,t+1);
end

% --- Constant ---
HDconst_big = zeros(nlag*nvar, nobs+1);
HDconst     = zeros(nvar, nobs+1);
if const > 0
    for t = 1:nobs
        k = kfun(t);
        CC = zeros(nlag*nvar, 1);
        CC(1:nvar) = F{k}(:,1);
        HDconst_big(:,t+1) = CC + Fcomp{k} * HDconst_big(:,t);
        HDconst(:,t+1) = Icomp * HDconst_big(:,t+1);
    end
end

% --- Linear trend ---
HDtrend_big = zeros(nlag*nvar, nobs+1);
HDtrend     = zeros(nvar, nobs+1);
if const > 1
    for t = 1:nobs
        k = kfun(t);
        TT = zeros(nlag*nvar, 1);
        TT(1:nvar) = F{k}(:,2);
        HDtrend_big(:,t+1) = TT*t + Fcomp{k} * HDtrend_big(:,t);
        HDtrend(:,t+1) = Icomp * HDtrend_big(:,t+1);
    end
end

% --- Quadratic trend ---
HDtrend2_big = zeros(nlag*nvar, nobs+1);
HDtrend2     = zeros(nvar, nobs+1);
if const > 2
    for t = 1:nobs
        k = kfun(t);
        TT2 = zeros(nlag*nvar, 1);
        TT2(1:nvar) = F{k}(:,3);
        HDtrend2_big(:,t+1) = TT2*(t^2) + Fcomp{k} * HDtrend2_big(:,t);
        HDtrend2(:,t+1) = Icomp * HDtrend2_big(:,t+1);
    end
end

% --- Sum: all components must add up to original data ---
HDendo = HDinit + HDconst + HDtrend + HDtrend2 + sum(HDshock,3);


%% Save and reshape all HDs
%==========================================================================
HD.shock = zeros(nobs+nlag, nvar, nvar);  % [nobs+nlag x shock x var]
for i = 1:nvar
    for j = 1:nvar
        HD.shock(:,j,i) = [nan(nlag,1); HDshock(i,2:end,j)'];
    end
end
HD.init   = [nan(nlag-1,nvar); HDinit(:,1:end)'];    % [nobs+nlag x nvar]
HD.const  = [nan(nlag,nvar);   HDconst(:,2:end)'];   % [nobs+nlag x nvar]
HD.trend  = [nan(nlag,nvar);   HDtrend(:,2:end)'];    % [nobs+nlag x nvar]
HD.trend2 = [nan(nlag,nvar);   HDtrend2(:,2:end)'];   % [nobs+nlag x nvar]
HD.endo   = [nan(nlag,nvar);   HDendo(:,2:end)'];     % [nobs+nlag x nvar]

% Update TVAR with structural impact matrices
for k = 1:nregimes
    TVAR.regime{k}.B = B{k};
end
