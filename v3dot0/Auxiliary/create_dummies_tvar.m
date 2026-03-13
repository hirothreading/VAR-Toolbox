function [yd, xd] = create_dummies_tvar(lambdaP, tauP, deltaP, epsilonP, nlag, muP, sigmaP, nvar)
%==========================================================================
% Create dummy observations implementing Minnesota-style priors for a
% Bayesian VAR. See Banbura, Giannone & Reichlin (2010, JAE).
%==========================================================================
% [yd, xd] = create_dummies_tvar(lambdaP, tauP, deltaP, epsilonP, nlag, muP, sigmaP, nvar)
% -------------------------------------------------------------------------
% INPUT
%   - lambdaP  : tightness on first lag (e.g. 0.2)
%   - tauP     : tightness on sum-of-coefficients (e.g. 10*lambdaP)
%   - deltaP   : (nvar x 1) prior mean of AR(1) coefficients
%   - epsilonP  : tightness on constant (e.g. 0.001)
%   - nlag     : number of lags
%   - muP      : (nvar x 1) sample mean of endogenous variables
%   - sigmaP   : (nvar x 1) residual std from univariate AR(1)
%   - nvar     : number of endogenous variables
% -------------------------------------------------------------------------
% OUTPUT
%   - yd : dummy dependent variable matrix
%   - xd : dummy regressor matrix
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit (create_dummies.m)
% =========================================================================

xd = [];
yd = [];
yd1 = [];
yd2 = [];
xd1 = [];
xd2 = [];

%% Dummy observations from equation (5) of Banbura et al. (2010)
if lambdaP > 0
    if epsilonP > 0
        yd1 = [diag(sigmaP .* deltaP) ./ lambdaP;
               zeros(nvar*(nlag-1), nvar);
               diag(sigmaP);
               zeros(1, nvar)];

        jp = diag(1:nlag);

        xd1 = [kron(jp, diag(sigmaP) ./ lambdaP),  zeros(nvar*nlag, 1);
               zeros(nvar, nvar*nlag+1);
               zeros(1, nvar*nlag),                  epsilonP];
    else
        yd1 = [diag(sigmaP .* deltaP) ./ lambdaP;
               zeros(nvar*(nlag-1), nvar);
               diag(sigmaP)];

        jp = diag(1:nlag);

        xd1 = [kron(jp, diag(sigmaP) ./ lambdaP);
               zeros(nvar, nvar*nlag)];
    end
end

%% Sum-of-coefficients dummy — equation (9) of Banbura et al. (2010)
if tauP > 0
    yd2 = diag(deltaP .* muP) ./ tauP;
    if epsilonP > 0
        xd2 = [kron(ones(1,nlag), yd2),  zeros(nvar, 1)];
    else
        xd2 = kron(ones(1,nlag), yd2);
    end
end

%% Concatenate
yd = [yd1; yd2];
xd = [xd1; xd2];
