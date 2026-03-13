function [beta, PROBLEM] = getcoef_tvar(mstar, sigma, ixx, maxtrys, nvar, nlag)
%==========================================================================
% Draw VAR coefficients from the conditional posterior, rejecting unstable
% draws. Used inside the Gibbs sampler of TVARmodel (Bayesian mode).
%==========================================================================
% [beta, PROBLEM] = getcoef_tvar(mstar, sigma, ixx, maxtrys, nvar, nlag)
% -------------------------------------------------------------------------
% INPUT
%   - mstar   : posterior mean of vec(beta), size (nvar*(nvar*nlag+1) x 1)
%   - sigma   : current covariance draw (nvar x nvar)
%   - ixx     : inv(X'X) matrix ((nvar*nlag+1) x (nvar*nlag+1))
%   - maxtrys : max number of draws before returning PROBLEM=1
%   - nvar    : number of endogenous variables
%   - nlag    : number of lags
% -------------------------------------------------------------------------
% OUTPUT
%   - beta    : draw from posterior (nvar*(nvar*nlag+1) x 1)
%   - PROBLEM : 1 if no stable draw found within maxtrys
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit (getcoef.m)
% =========================================================================

PROBLEM = 0;
ncoeff = nvar * nlag + 1;  % coefficients per equation (lags + constant)
vstar = kron(sigma, ixx);

check = -1;
tryx = 1;
while check < 0 && tryx < maxtrys
    beta = mstar + (randn(1, nvar*ncoeff) * chol(vstar))';

    % Check stability: eigenvalues of companion matrix must be < 1
    S = stability_tvar(beta, nvar, nlag);
    if S == 0
        check = 10;  % stable draw found
    else
        tryx = tryx + 1;
    end
end

if S > 0
    PROBLEM = 1;
end
