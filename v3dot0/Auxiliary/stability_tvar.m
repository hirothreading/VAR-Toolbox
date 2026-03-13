function S = stability_tvar(beta, nvar, nlag)
%==========================================================================
% Check VAR stability by testing whether the maximum eigenvalue of the
% companion matrix is less than 1.
%==========================================================================
% S = stability_tvar(beta, nvar, nlag)
% -------------------------------------------------------------------------
% INPUT
%   - beta : vectorised coefficient vector (nvar*(nvar*nlag+1) x 1)
%   - nvar : number of endogenous variables
%   - nlag : number of lags
% -------------------------------------------------------------------------
% OUTPUT
%   - S : 0 if stable (max|eigenvalue| < 1), 1 if unstable
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit (stability.m)
% =========================================================================

ncoeff = nvar * nlag + 1;  % lags + constant per equation

% Build companion matrix
FF = zeros(nvar*nlag, nvar*nlag);
if nlag > 1
    FF(nvar+1:nvar*nlag, 1:nvar*(nlag-1)) = eye(nvar*(nlag-1));
end

% Extract lag coefficients (exclude constant) and place in top block
temp = reshape(beta, ncoeff, nvar);
temp = temp(1:nvar*nlag, :)';        % nvar x (nvar*nlag)
FF(1:nvar, 1:nvar*nlag) = temp;

% Stability check
ee = max(abs(eig(FF)));
S = (ee >= 1);
