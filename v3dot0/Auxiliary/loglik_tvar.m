function out = loglik_tvar(beta_mat, sigma, Y, X)
%==========================================================================
% Evaluate the Gaussian log-likelihood for a VAR with given coefficients
% and covariance matrix.
%==========================================================================
% out = loglik_tvar(beta_mat, sigma, Y, X)
% -------------------------------------------------------------------------
% INPUT
%   - beta_mat : coefficient matrix ((nvar*nlag+1) x nvar)
%   - sigma    : covariance matrix (nvar x nvar)
%   - Y        : dependent variable (T x nvar)
%   - X        : regressors (T x (nvar*nlag+1))
% -------------------------------------------------------------------------
% OUTPUT
%   - out : log-likelihood value (scalar)
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit (loglik.m)
% =========================================================================

T = size(Y, 1);
V = Y - X * beta_mat;   % residuals (T x nvar)

% Use Cholesky for numerically stable inverse and log-determinant
[L_chol, flag] = chol(sigma, 'lower');
if flag ~= 0
    out = -inf;
    return;
end
isigma = L_chol' \ (L_chol \ eye(size(sigma,1)));
dsigma = 2 * sum(log(diag(L_chol)));  % log|sigma|

% Sum of quadratic forms: sum_t v_t' * inv(sigma) * v_t
sterm = sum(sum((V * isigma) .* V));   % vectorised, avoids loop

out = (T/2) * (-dsigma) - 0.5 * sterm;
