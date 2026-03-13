function [VD, TVAR] = TVARvd(TVAR, VARopt)
%==========================================================================
% Compute regime-specific forecast error variance decompositions for a
% TVAR estimated with TVARmodel. Delegates to VARvd per regime.
%
% Not available for 'iv' identification.
%==========================================================================
% [VD, TVAR] = TVARvd(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel.m
%   - VARopt : options (ident, nsteps, etc.)
% -------------------------------------------------------------------------
% OUTPUT
%   - VD     : (nsteps x nshocks x nvar x nregimes) variance decomposition
%   - TVAR   : updated (regime{k}.B populated)
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------

if strcmp(VARopt.ident, 'iv')
    error('TVARvd: FEVD not available with external instruments (iv) identification');
end

nvar     = TVAR.nvar;
nsteps   = VARopt.nsteps;
nregimes = TVAR.nregimes;

VD = nan(nsteps, nvar, nvar, nregimes);

for k = 1:nregimes
    [VD_k, regime_k] = VARvd(TVAR.regime{k}, VARopt);
    VD(:,:,:,k) = VD_k;
    TVAR.regime{k} = regime_k;
end
