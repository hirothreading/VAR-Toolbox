function TVARprint(TVAR, TVARopt, approx)
%==========================================================================
% Print the output of a TVAR estimation to screen
%==========================================================================
% TVARprint(TVAR, TVARopt, approx)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR    : structure, result of TVARmodel.m
%   - TVARopt : options of the TVAR (result of TVARmodel.m)
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - approx  : number of decimal digits  [default = 4]
% -------------------------------------------------------------------------
% EXAMPLE
%   - See TVARToolbox_Primer.m in "../Primer/"
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================


%% Check inputs
%--------------------------------------------------------------------------
if ~exist('TVAR','var')
    error('You need to provide TVAR structure, result of TVARmodel');
end
if ~exist('TVARopt','var')
    error('You need to provide TVARopt structure');
end
if ~exist('approx','var') || isempty(approx)
    approx = 4;
end

vnames   = TVARopt.vnames;
rnames   = TVARopt.rnames;
nvar     = TVAR.nvar;
nlag     = TVAR.nlag;
const    = TVAR.const;

if isempty(vnames)
    % Fall back to generic names if not set
    vnames = {};
    for ii = 1:nvar
        vnames{ii} = ['y' num2str(ii)];
    end
end


%% Header
%--------------------------------------------------------------------------
disp(' ')
disp('=================================================================')
disp('  Threshold VAR (TVAR) Estimation Results')
disp('=================================================================')
fprintf('  Threshold variable : ');
if TVAR.thrvar_idx == 0
    fprintf('external\n');
else
    if ~isempty(vnames) && length(vnames) >= TVAR.thrvar_idx
        fprintf('%s (column %d)\n', vnames{TVAR.thrvar_idx}, TVAR.thrvar_idx);
    else
        fprintf('column %d\n', TVAR.thrvar_idx);
    end
end
fprintf('  Delay              : %d\n', TVAR.delay);
fprintf('  Estimated threshold: %.*f\n', approx, TVAR.thresh);
fprintf('  Total observations : %d\n', TVAR.nobs);
fprintf('  Lag order          : %d\n', TVAR.nlag);
switch const
    case 0; fprintf('  Deterministics     : none\n');
    case 1; fprintf('  Deterministics     : constant\n');
    case 2; fprintf('  Deterministics     : constant + trend\n');
    case 3; fprintf('  Deterministics     : constant + trend + trend^2\n');
end
disp(' ')


%% Regime summary table
%--------------------------------------------------------------------------
disp('-----------------------------------------------------------------')
disp('  Regime summary')
disp('-----------------------------------------------------------------')
fprintf('  %-15s  %6s  %8s  %10s  %8s\n', 'Regime', 'Obs', '%Sample', 'Max |Eig|', 'Stable?');
fprintf('  %-15s  %6s  %8s  %10s  %8s\n', '------', '---', '-------', '---------', '-------');
for k = 1:TVAR.nregimes
    rname = rnames{k};
    nk    = TVAR.regime{k}.nobs;
    pct   = 100 * nk / TVAR.nobs;
    eig_k = TVAR.regime{k}.maxEig;
    stab  = 'Yes';
    if eig_k >= 1; stab = 'NO'; end
    fprintf('  %-15s  %6d  %7.1f%%  %10.4f  %8s\n', rname, nk, pct, eig_k, stab);
end
disp(' ')


%% Build row labels for coefficient table (same logic as VARprint)
%--------------------------------------------------------------------------
switch const
    case 0; det_labels = {};
    case 1; det_labels = {'c'};
    case 2; det_labels = {'c'; 'trend'};
    case 3; det_labels = {'c'; 'trend'; 'trend2'};
end

row_labels = det_labels(:);
for jj = 1:nlag
    for ii = 1:nvar
        row_labels{end+1} = [vnames{ii} '(-' num2str(jj) ')']; %#ok<AGROW>
    end
end

% Exogenous labels
if TVAR.regime{1}.nvar_ex > 0 && ~isempty(TVARopt.vnames_ex)
    vnames_ex = TVARopt.vnames_ex;
    for ii = 1:TVAR.regime{1}.nvar_ex
        row_labels{end+1} = vnames_ex{ii}; %#ok<AGROW>
    end
    for jj = 1:TVAR.regime{1}.nlag_ex
        for ii = 1:TVAR.regime{1}.nvar_ex
            row_labels{end+1} = [vnames_ex{ii} '(-' num2str(jj) ')']; %#ok<AGROW>
        end
    end
end


%% Print regime-specific coefficients
%--------------------------------------------------------------------------
for k = 1:TVAR.nregimes
    disp('-----------------------------------------------------------------')
    fprintf('  %s  |  Coefficient Matrix (Ft)\n', rnames{k});
    disp('-----------------------------------------------------------------')

    info.cnames = char(vnames);
    info.rnames = char([{' '}; row_labels]);
    mprint(TVAR.regime{k}.Ft, info);

    disp(' ')
    fprintf('  %s  |  Residual Covariance Matrix\n', rnames{k});
    disp(TVAR.regime{k}.sigma)
end


%% R-squared by regime
%--------------------------------------------------------------------------
disp('-----------------------------------------------------------------')
disp('  R-squared by regime and equation')
disp('-----------------------------------------------------------------')
fprintf('  %-15s', '');
for ii = 1:nvar
    fprintf('  %10s', vnames{ii});
end
fprintf('\n');
for k = 1:TVAR.nregimes
    fprintf('  %-15s', rnames{k});
    for ii = 1:nvar
        eval(['rsqr = TVAR.regime{k}.eq' num2str(ii) '.rsqr;']);
        fprintf('  %10.4f', rsqr);
    end
    fprintf('\n');
end
disp(' ')
