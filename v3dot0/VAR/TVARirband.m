function [INF, SUP, MED, BAR] = TVARirband(TVAR, VARopt)
%==========================================================================
% Compute bootstrap confidence bands for regime-specific impulse response
% functions from a Threshold VAR estimated with TVARmodel / TVARir.
%
% The bootstrap uses a fixed-design approach (standard in the TVAR
% literature): artificial data is generated with the same threshold
% variable values as observed, so the regime assignments at each time step
% are fixed at their observed values when propagating shocks. The threshold
% gamma is re-estimated in each bootstrap draw.
%==========================================================================
% [INF, SUP, MED, BAR] = TVARirband(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure, result of TVARmodel.m (and optionally TVARir.m
%              for 'sign' identification where regime{k}.B is pre-set)
%   - VARopt : VAR options.  Key fields used:
%                .nsteps   - IRF horizon
%                .ndraws   - number of bootstrap draws
%                .pctg     - confidence level  (e.g. 68 or 95)
%                .method   - 'bs' (standard) or 'wild'
%                .ident    - identification scheme (same as TVARir)
%                .mult     - print-every-N-draws counter
% -------------------------------------------------------------------------
% OUTPUT
%   - INF(h,n,s,k) : lower confidence band
%   - SUP(h,n,s,k) : upper confidence band
%   - MED(h,n,s,k) : median IRF
%   - BAR(h,n,s,k) : mean IRF
%   Dimensions: (nsteps x nvar x nvar x nregimes)
%   k=1 => regime 1,  k=2 => regime 2
% -------------------------------------------------------------------------
% NOTE
%   For 'iv' identification TVAR.IV must be set.
%   For 'sign' identification the B matrix for each regime must be set
%   via SR before calling TVARirband (sign restrictions are NOT
%   re-identified in each bootstrap draw in this implementation).
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
if ~exist('VARopt','var')
    error('You need to provide VARopt structure');
end

nsteps   = VARopt.nsteps;
ndraws   = VARopt.ndraws;
pctg     = VARopt.pctg;
method   = VARopt.method;
ident    = VARopt.ident;
nvar     = TVAR.nvar;
nlag     = TVAR.nlag;
const    = TVAR.const;
nobse    = TVAR.nobs;
nobs     = size(TVAR.ENDO, 1);
nregimes = TVAR.nregimes;
ENDO     = TVAR.ENDO;
EXOG     = TVAR.EXOG;
regime_idx  = TVAR.regime_idx;
thrvar_idx  = TVAR.thrvar_idx;
delay       = TVAR.delay;

% Regime residuals and coefficients
resid1 = TVAR.regime{1}.resid;   % (nobs1 x nvar)
resid2 = TVAR.regime{2}.resid;   % (nobs2 x nvar)
nobs1  = TVAR.regime{1}.nobs;
nobs2  = TVAR.regime{2}.nobs;
Ft1    = TVAR.regime{1}.Ft;
Ft2    = TVAR.regime{2}.Ft;

% Check iv instrument
if strcmp(ident,'iv') && isempty(TVAR.IV)
    error('TVARirband: VARopt.ident=''iv'' requires TVAR.IV to be set');
end

% Pre-allocate IRF storage: [nsteps x nvar x nvar x nregimes x ndraws]
IR_store = nan(nsteps, nvar, nvar, nregimes, ndraws);
MED      = zeros(nsteps, nvar, nvar, nregimes);
BAR      = zeros(nsteps, nvar, nvar, nregimes);

% Artificial data buffer
y_art = zeros(nobs, nvar);


%% Bootstrap loop
%--------------------------------------------------------------------------
tt = 1;   % accepted draws
ww = 1;   % display counter

while tt <= ndraws

    % Progress display
    if tt == VARopt.mult * ww
        disp(['TVARirband: ' num2str(tt) ' / ' num2str(ndraws) ' draws']);
        ww = ww + 1;
    end

    %% STEP 1: draw residuals within each regime
    if strcmp(method, 'bs')
        u1 = resid1(ceil(nobs1 * rand(nobs1,1)), :);
        u2 = resid2(ceil(nobs2 * rand(nobs2,1)), :);

    elseif strcmp(method, 'wild')
        rr1 = 1 - 2*(rand(nobs1,1) > 0.5);
        rr2 = 1 - 2*(rand(nobs2,1) > 0.5);
        u1  = resid1 .* (rr1 * ones(1,nvar));
        u2  = resid2 .* (rr2 * ones(1,nvar));

        % Wild bootstrap for instrument: apply same perturbation to IV rows
        if strcmp(ident,'iv')
            IV_full = TVAR.IV;
            rr_full = nan(nobse,1);
            rr_full(regime_idx==1) = rr1;
            rr_full(regime_idx==2) = rr2;
            % IV is indexed from row nlag+1:nobs in original space
            rr_iv = rr_full;   % same length as Y rows
            IV_boot = IV_full;
            IV_boot(nlag+1:end, :) = IV_full(nlag+1:end,:) .* ...
                (rr_iv * ones(1,size(IV_full,2)));
        end
    else
        error(['TVARirband: unknown bootstrap method ''' method '''']);
    end

    % Reassemble into full-sample order (fixed-design: regime assignments unchanged)
    u = zeros(nobse, nvar);
    u(regime_idx==1, :) = u1;
    u(regime_idx==2, :) = u2;


    %% STEP 2: generate artificial data
    % Seed first nlag rows from actual data
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
    % Append exogenous if present (same as VARirband)
    if ~isempty(EXOG) && TVAR.regime{1}.nvar_ex > 0
        LAGplus = [LAGplus TVAR.X_EX(1,:)];
    end

    % Fixed-design: use Ft_k selected by the original regime_idx
    for jj = nlag+1 : nobs
        t_aligned = jj - nlag;  % row index in Y (1-based)
        k_t = regime_idx(t_aligned);
        if k_t == 1; Ft_t = Ft1; else; Ft_t = Ft2; end
        for mm = 1:nvar
            y_art(jj,mm) = LAGplus * Ft_t(:,mm) + u(t_aligned, mm);
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
            if ~isempty(EXOG) && TVAR.regime{1}.nvar_ex > 0
                LAGplus = [LAGplus TVAR.X_EX(tidx,:)]; %#ok<AGROW>
            end
        end
    end


    %% STEP 3: re-estimate TVAR on artificial data
    try
        if ~isempty(EXOG)
            TVAR_draw = TVARmodel(y_art, nlag, const, thrvar_idx, delay, [], EXOG, TVAR.regime{1}.nlag_ex);
        else
            TVAR_draw = TVARmodel(y_art, nlag, const, thrvar_idx, delay);
        end
    catch
        % Threshold estimation failed on this draw — skip
        continue
    end

    % Both regimes must be stable
    if TVAR_draw.regime{1}.maxEig >= 0.9999 || TVAR_draw.regime{2}.maxEig >= 0.9999
        continue
    end

    % Attach instrument for iv identification
    if strcmp(ident,'iv')
        if strcmp(method,'wild')
            TVAR_draw.IV = IV_boot;
        else
            TVAR_draw.IV = TVAR.IV;
        end
    end


    %% STEP 4: compute IRFs on this draw
    try
        [IR_draw, ~] = TVARir(TVAR_draw, VARopt);
    catch
        continue
    end

    IR_store(:,:,:,:,tt) = IR_draw;
    tt = tt + 1;
end

disp('TVARirband: -- Done!');
disp(' ');


%% Compute confidence bands
%--------------------------------------------------------------------------
npctg = numel(pctg);
if npctg > 1
    INF = zeros(nsteps, nvar, nvar, nregimes, npctg);
    SUP = zeros(nsteps, nvar, nvar, nregimes, npctg);
else
    INF = zeros(nsteps, nvar, nvar, nregimes);
    SUP = zeros(nsteps, nvar, nvar, nregimes);
end

for k = 1:nregimes
    if npctg > 1
        for pp = 1:npctg
            INF(:,:,:,k,pp) = prctile(IR_store(:,:,:,k,:), (100-pctg(pp))/2, 5);
            SUP(:,:,:,k,pp) = prctile(IR_store(:,:,:,k,:), 100-(100-pctg(pp))/2, 5);
        end
    else
        INF(:,:,:,k) = prctile(IR_store(:,:,:,k,:), (100-pctg)/2, 5);
        SUP(:,:,:,k) = prctile(IR_store(:,:,:,k,:), 100-(100-pctg)/2, 5);
    end
    MED(:,:,:,k) = prctile(IR_store(:,:,:,k,:), 50, 5);
    BAR(:,:,:,k) = mean(IR_store(:,:,:,k,:), 5);
end
