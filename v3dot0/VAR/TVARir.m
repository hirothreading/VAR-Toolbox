function [IR, TVAR] = TVARir(TVAR, VARopt)
%==========================================================================
% Compute Generalized Impulse Response Functions (GIRFs) for a Threshold
% VAR estimated with TVARmodel, following Koop, Pesaran & Potter (1996).
%
% GIRFs are computed by forward simulation, allowing the system to switch
% regimes endogenously during the IRF horizon. The same random innovations
% are used for shocked and unshocked paths (variance reduction).
%
% Identification:
%   'short' — Cholesky decomposition of regime-specific covariance
%
% Output IR is 4-D: (nsteps x nvar x nshocks x nregimes).
%==========================================================================
% [IR, TVAR] = TVARir(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel.m
%   - VARopt : VAR options (ident, nsteps, nreps_girf, shock_scale, etc.)
% -------------------------------------------------------------------------
% OUTPUT
%   - IR     : (nsteps x nvar x nshocks x nregimes) impulse responses
%   - TVAR   : updated structure
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('TVAR','var')
    error('You need to provide TVAR structure, result of TVARmodel');
end
if ~exist('VARopt','var')
    error('You need to provide VARopt structure');
end

ident      = VARopt.ident;
nvar       = TVAR.nvar;
nlag       = TVAR.nlag;
nsteps     = VARopt.nsteps;
nregimes   = TVAR.nregimes;

% GIRF-specific options
if isfield(VARopt, 'nreps_girf')
    nreps_sim = VARopt.nreps_girf;
else
    nreps_sim = 50;
end
if isfield(VARopt, 'shock_scale')
    shock_scale = VARopt.shock_scale;
else
    shock_scale = 1;
end

% Threshold specification
delay        = TVAR.delay;
thrvar_idx   = VARopt.thrvar_idx;
tartransform = VARopt.tartransform;


%% Get coefficients in Mumtaz ordering [lags|constant]
%==========================================================================
% For point estimate GIRFs, convert Ft from VARmodel ordering
% [const|lags] to Mumtaz ordering [lags|const]
for k = 1:nregimes
    Ftk = TVAR.regime{k}.Ft;  % (ntotcoeff x nvar)
    % VARmodel: row 1 = constant, rows 2:end = lags
    % Mumtaz: rows 1:nvar*nlag = lags, last row = constant
    beta_mumtaz = [Ftk(2:nvar*nlag+1, :); Ftk(1, :)];
    TVAR.regime{k}.beta_mumtaz = beta_mumtaz(:);  % vectorised
end

beta1 = TVAR.regime{1}.beta_mumtaz;
beta2 = TVAR.regime{2}.beta_mumtaz;
sigma1 = TVAR.regime{1}.sigma;
sigma2 = TVAR.regime{2}.sigma;
thresh = TVAR.thresh;


%% Compute GIRFs
%==========================================================================
if strcmp(ident, 'short')
    % Cholesky identification via GIRF simulation
    shock_pos = 1:nvar;  % all shocks

    GIRF = TVARgirf(beta1, sigma1, beta2, sigma2, thresh, ...
        TVAR.Y, TVAR.Ystar, nvar, nlag, delay, nsteps, nreps_sim, ...
        shock_pos, shock_scale, tartransform, thrvar_idx);

    IR = nan(nsteps, nvar, nvar, nregimes);
    IR(:,:,:,1) = GIRF.regime1;
    IR(:,:,:,2) = GIRF.regime2;

else
    error('TVARir: ident=''%s'' is not yet supported for TVAR GIRFs. Use ''short''.', ident);
end


%% Store identification info
%==========================================================================
for k = 1:nregimes
    TVAR.regime{k}.B = chol(TVAR.regime{k}.sigma)';  % lower triangular
end

end
