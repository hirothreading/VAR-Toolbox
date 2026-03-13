function TVARprint(TVAR, VARopt, approx)
%==========================================================================
% Print the estimation results of a Threshold VAR.
%==========================================================================
% TVARprint(TVAR, VARopt, approx)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel.m
%   - VARopt : options with vnames populated
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - approx : decimal digits for display [default = 4]
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------

if ~exist('approx','var')
    approx = 4;
end

vnames = VARopt.vnames;
if isempty(vnames)
    error('You need to add labels for endogenous variables in VARopt.vnames');
end

nvar = TVAR.nvar;
nlag = TVAR.nlag;

%% Header
%==========================================================================
disp(' ');
disp('============================================================');
disp('  THRESHOLD VAR ESTIMATION RESULTS');
disp('============================================================');
disp(['  Method:       ' upper(TVAR.method)]);
disp(['  Variables:    ' num2str(nvar)]);
disp(['  Lags:         ' num2str(nlag)]);
disp(['  Observations: ' num2str(TVAR.nobs)]);
disp(['  Threshold:    ' num2str(TVAR.thresh, approx+2)]);
disp(['  Delay:        ' num2str(TVAR.delay)]);
disp(' ');

%% Bayesian diagnostics
%==========================================================================
if strcmp(TVAR.method, 'bayes')
    disp('  --- MCMC Diagnostics ---');
    disp(['  Total replications:  ' num2str(TVAR.nreps)]);
    disp(['  Burn-in:             ' num2str(TVAR.nburn)]);
    disp(['  Thinning:            ' num2str(TVAR.nskip)]);
    disp(['  Retained draws:      ' num2str(TVAR.ndraws)]);
    disp(['  Acceptance rate:     ' num2str(TVAR.acceptance_rate, '%.3f')]);
    disp(['  Threshold std:       ' num2str(std(TVAR.thresh_draws), '%.4f')]);
    disp(' ');
end

%% Regime results
%==========================================================================
for k = 1:TVAR.nregimes
    reg = TVAR.regime{k};
    if isfield(VARopt, 'regime_names') && ~isempty(VARopt.regime_names)
        rname = VARopt.regime_names{k};
    else
        rname = ['Regime ' num2str(k)];
    end

    disp('------------------------------------------------------------');
    disp(['  ' rname ' (' num2str(reg.nobs) ' observations, ' ...
        num2str(100*reg.nobs/TVAR.nobs, '%.1f') '% of sample)']);
    disp('------------------------------------------------------------');

    % Coefficient matrix
    info.cnames = char(vnames);
    % Build row labels
    vtext = {};
    if TVAR.const >= 1; vtext{end+1} = 'c'; end
    if TVAR.const >= 2; vtext{end+1} = 'trend'; end
    if TVAR.const >= 3; vtext{end+1} = 'trend2'; end
    for jj = 1:nlag
        for ii = 1:nvar
            vtext{end+1} = [vnames{ii} '(-' num2str(jj) ')'];
        end
    end
    info.rnames = char([{' '}; vtext']);

    disp(' ');
    disp('  Coefficients:');
    mprint(reg.Ft, info);

    disp(' ');
    disp('  Eigenvalues:');
    disp(eig(reg.Fcomp));

    disp('  Covariance matrix:');
    disp(reg.sigma);

    % R-squared per equation (if available)
    if isfield(reg, 'eq')
        fprintf('  R-squared: ');
        for j = 1:nvar
            fprintf('%s=%.3f  ', vnames{j}, reg.eq{j}.rsqr);
        end
        fprintf('\n');
    end
    disp(' ');
end

disp('============================================================');
