function out = TVARtest(ENDO, nlag, const, thrvar_idx, delay, VARopt, THRVAR_EX, EXOG, nlag_ex)
%==========================================================================
% Bootstrap test of linearity (no threshold) against a two-regime
% Threshold VAR, following Hansen (1999).
%
% H0: linear VAR (single regime)
% H1: two-regime TVAR with one threshold
%
% The test statistic is an F-type statistic:
%   F = n * (RSS_linear - RSS_tvar) / RSS_tvar
%
% Under H0 the threshold gamma is not identified, so standard asymptotic
% distributions are non-standard. Bootstrap p-values are computed by
% generating data under H0 (linear VAR) and re-running the full grid
% search on each bootstrap draw.
%==========================================================================
% out = TVARtest(ENDO, nlag, const, thrvar_idx, delay, VARopt, THRVAR_EX, EXOG, nlag_ex)
% -------------------------------------------------------------------------
% INPUT
%   - ENDO       : (nobs x nvar) matrix of endogenous variables
%   - nlag       : lag order
%   - const      : 0=none; 1=constant; 2=const+trend; 3=const+trend+trend^2
%   - thrvar_idx : column index in ENDO for threshold variable (0=external)
%   - delay      : delay d for threshold variable  [default = 1]
%   - VARopt     : VAR options structure (uses VARopt.ndraws, VARopt.method,
%                  VARopt.mult for display)
% -------------------------------------------------------------------------
% OPTIONAL INPUT
%   - THRVAR_EX  : (nobs x 1) external threshold variable (if thrvar_idx=0)
%   - EXOG       : (nobs x nvar_ex) exogenous variables
%   - nlag_ex    : lag order for exogenous variables
% -------------------------------------------------------------------------
% OUTPUT
%   - out.F_stat     : observed F-type test statistic
%   - out.pval       : bootstrap p-value
%   - out.cv         : [cv90, cv95, cv99] bootstrap critical values
%   - out.nboot      : number of accepted bootstrap draws
%   - out.RSS_linear : RSS of the linear VAR
%   - out.RSS_tvar   : RSS of the TVAR at the optimal threshold
%   - out.grid_gamma : threshold grid from TVARmodel
%   - out.grid_RSS   : TVAR RSS at each grid point
% -------------------------------------------------------------------------
% EXAMPLE
%   - See TVARToolbox_Primer.m in "../Primer/"
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
% Reference: Hansen, B.E. (1999). "Testing for Linearity." Journal of
%            Economic Surveys, 13(5), 551-576.
%==========================================================================


%% Check inputs
%--------------------------------------------------------------------------
if ~exist('delay','var') || isempty(delay)
    delay = 1;
end
if ~exist('const','var') || isempty(const)
    const = 1;
end
if ~exist('VARopt','var')
    VARopt = VARoption;
end

nboot  = VARopt.ndraws;
method = VARopt.method;

has_threx  = exist('THRVAR_EX','var') && ~isempty(THRVAR_EX);
has_exog   = exist('EXOG','var') && ~isempty(EXOG);
if has_exog && (~exist('nlag_ex','var') || isempty(nlag_ex))
    nlag_ex = 0;
end

[nobs, nvar] = size(ENDO);


%% Step 1: estimate linear VAR and TVAR on actual data
%--------------------------------------------------------------------------
if has_exog
    [VAR0, ~] = VARmodel(ENDO, nlag, const, EXOG, nlag_ex);
else
    [VAR0, ~] = VARmodel(ENDO, nlag, const);
end

% RSS of linear VAR
resid0  = VAR0.resid;
RSS0    = sum(sum(resid0.^2));

% TVAR: run TVARmodel
if thrvar_idx == 0 && has_threx
    if has_exog
        [TVAR, ~] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, THRVAR_EX, EXOG, nlag_ex);
    else
        [TVAR, ~] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, THRVAR_EX);
    end
else
    if has_exog
        [TVAR, ~] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay, [], EXOG, nlag_ex);
    else
        [TVAR, ~] = TVARmodel(ENDO, nlag, const, thrvar_idx, delay);
    end
end

RSS1     = TVAR.RSS;
nobse    = TVAR.nobs;
F_stat   = nobse * (RSS0 - RSS1) / RSS1;


%% Step 2: bootstrap under H0 (linear VAR)
%--------------------------------------------------------------------------
Ft0      = VAR0.Ft;        % (ntotcoeff x nvar) linear VAR coefficients
resid0   = VAR0.resid;     % (nobse x nvar) linear VAR residuals

% Grid options for the bootstrap TVAR (re-use TVARopt from TVAR)
thrvar_idx_boot = thrvar_idx;
delay_boot      = delay;

F_b  = nan(nboot, 1);
tt   = 1;   % accepted draws
ww   = 1;   % display counter

y_art = zeros(nobs, nvar);

while tt <= nboot

    % Progress display
    if tt == VARopt.mult * ww
        disp(['TVARtest bootstrap: ' num2str(tt) ' / ' num2str(nboot)]);
        ww = ww + 1;
    end

    %% Draw residuals under H0
    if strcmp(method, 'bs')
        u = resid0(ceil(nobse * rand(nobse,1)), :);
    elseif strcmp(method, 'wild')
        rr = 1 - 2*(rand(nobse,1) > 0.5);
        u  = resid0 .* (rr * ones(1,nvar));
    else
        error(['TVARtest: unknown bootstrap method ''' method '''']);
    end

    %% Generate artificial data under H0 (single regime)
    % Initialise with actual data for the first nlag rows
    LAG = [];
    for jj = 1:nlag
        y_art(jj,:) = ENDO(jj,:);
        LAG = [y_art(jj,:) LAG]; %#ok<AGROW>
    end
    T = (1:nobse)';
    switch const
        case 0; LAGplus = LAG;
        case 1; LAGplus = [1 LAG];
        case 2; LAGplus = [1 T(1) LAG];
        case 3; LAGplus = [1 T(1) T(1).^2 LAG];
    end

    for jj = nlag+1 : nobs
        for mm = 1:nvar
            y_art(jj,mm) = LAGplus * Ft0(:,mm) + u(jj-nlag, mm);
        end
        if jj < nobs
            LAG = [y_art(jj,:) LAG(1, 1:(nlag-1)*nvar)]; %#ok<AGROW>
            tidx = jj - nlag + 1;
            switch const
                case 0; LAGplus = LAG;
                case 1; LAGplus = [1 LAG];
                case 2; LAGplus = [1 T(tidx) LAG];
                case 3; LAGplus = [1 T(tidx) T(tidx).^2 LAG];
            end
        end
    end

    %% Re-estimate linear VAR and TVAR on artificial data
    if has_exog
        VAR0_b = VARmodel(y_art, nlag, const, EXOG, nlag_ex);
    else
        VAR0_b = VARmodel(y_art, nlag, const);
    end
    RSS0_b = sum(sum(VAR0_b.resid.^2));

    try
        if thrvar_idx == 0 && has_threx
            if has_exog
                TVAR_b = TVARmodel(y_art, nlag, const, thrvar_idx_boot, delay_boot, THRVAR_EX, EXOG, nlag_ex);
            else
                TVAR_b = TVARmodel(y_art, nlag, const, thrvar_idx_boot, delay_boot, THRVAR_EX);
            end
        else
            if has_exog
                TVAR_b = TVARmodel(y_art, nlag, const, thrvar_idx_boot, delay_boot, [], EXOG, nlag_ex);
            else
                TVAR_b = TVARmodel(y_art, nlag, const, thrvar_idx_boot, delay_boot);
            end
        end
        RSS1_b = TVAR_b.RSS;
    catch
        % If threshold estimation fails (e.g., degenerate sample), skip
        continue
    end

    nobse_b   = TVAR_b.nobs;
    F_b(tt)   = nobse_b * (RSS0_b - RSS1_b) / RSS1_b;
    tt = tt + 1;
end

disp('TVARtest bootstrap: Done!');
disp(' ');

% Remove any NaN entries (shouldn't happen, but guard)
F_b = F_b(~isnan(F_b));


%% Compute p-value and critical values
%--------------------------------------------------------------------------
pval  = mean(F_b >= F_stat);
cv90  = quantile(F_b, 0.90);
cv95  = quantile(F_b, 0.95);
cv99  = quantile(F_b, 0.99);


%% Display results
%--------------------------------------------------------------------------
disp('=================================================================')
disp('  TVAR Linearity Test (Hansen 1999)')
disp('=================================================================')
fprintf('  H0: linear VAR   vs.   H1: two-regime TVAR\n\n');
fprintf('  F-statistic      : %8.4f\n', F_stat);
fprintf('  Bootstrap p-value: %8.4f  (%d draws)\n', pval, length(F_b));
fprintf('  Critical values  :  90%% = %6.3f   95%% = %6.3f   99%% = %6.3f\n', cv90, cv95, cv99);
fprintf('  RSS (linear VAR) : %8.4f\n', RSS0);
fprintf('  RSS (TVAR)       : %8.4f\n', RSS1);
fprintf('  Threshold        : %8.4f\n', TVAR.thresh);
disp('=================================================================')
disp(' ')


%% Assemble output struct
%--------------------------------------------------------------------------
out.F_stat      = F_stat;
out.pval        = pval;
out.cv          = [cv90, cv95, cv99];
out.nboot       = length(F_b);
out.RSS_linear  = RSS0;
out.RSS_tvar    = RSS1;
out.thresh      = TVAR.thresh;
out.grid_gamma  = TVAR.grid_gamma;
out.grid_RSS    = TVAR.grid_RSS;
out.F_boot      = F_b;   % full bootstrap distribution (useful for plotting)
