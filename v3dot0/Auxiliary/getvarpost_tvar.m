function [post, lik1, lik2, prior] = getvarpost_tvar(Y, X, beta1, beta2, sigma1, sigma2, nlag, tar, tarmean, tarvariance, Ystar, ncrit)
%==========================================================================
% Evaluate the log-posterior of the threshold value for the TVAR model.
%   log p(tau | data) = log L(regime 1) + log L(regime 2) + log p(tau)
%==========================================================================
% [post, lik1, lik2, prior] = getvarpost_tvar(...)
% -------------------------------------------------------------------------
% INPUT
%   - Y, X       : full-sample data (T x nvar), (T x (nvar*nlag+1))
%   - beta1,beta2: vectorised coefficients per regime (nvar*(nvar*nlag+1) x 1)
%   - sigma1,sigma2: covariance matrices per regime (nvar x nvar)
%   - nlag       : lag order
%   - tar        : candidate threshold value (scalar)
%   - tarmean    : prior mean of threshold
%   - tarvariance: prior variance of threshold
%   - Ystar      : threshold variable vector (T x 1)
%   - ncrit      : minimum number of observations per regime
% -------------------------------------------------------------------------
% OUTPUT
%   - post  : log-posterior (scalar)
%   - lik1  : log-likelihood in regime 1
%   - lik2  : log-likelihood in regime 2
%   - prior : log-prior on threshold
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit (getvarpost.m)
% =========================================================================

nvar = size(Y, 2);
ncoeff = nvar * nlag + 1;

% Split sample by threshold
e1 = (Ystar <= tar);
e2 = (Ystar > tar);

% Check minimum observations per regime
if sum(e1) < ncrit || sum(e2) < ncrit
    post = -inf;
    lik1 = -inf;
    lik2 = -inf;
    prior = -inf;
    return;
end

Y1 = Y(e1,:);  X1 = X(e1,:);
Y2 = Y(e2,:);  X2 = X(e2,:);

% Regime log-likelihoods
beta1_mat = reshape(beta1, ncoeff, nvar);
beta2_mat = reshape(beta2, ncoeff, nvar);
lik1 = loglik_tvar(beta1_mat, sigma1, Y1, X1);
lik2 = loglik_tvar(beta2_mat, sigma2, Y2, X2);

% Log-prior on threshold: N(tarmean, tarvariance)
res = tar - tarmean;
prior = -0.5 * log(2*pi*tarvariance) - 0.5 * (res^2) / tarvariance;

% Log-posterior
post = lik1 + lik2 + prior;
if isinf(post) || ~isreal(post) || isnan(post)
    post = -inf;
end
