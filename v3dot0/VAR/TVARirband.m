function [INF, SUP, MED, BAR] = TVARirband(TVAR, VARopt)
%==========================================================================
% Compute confidence / credible bands for TVAR Generalized IRFs.
%
%   Frequentist (tvar_method='freq'): bootstrap with residual resampling,
%     GIRF re-computed at each bootstrap draw.
%   Bayesian   (tvar_method='bayes'): GIRF computed at each posterior draw.
%
% Output arrays are 4-D (or 5-D with multiple pctg):
%   (nsteps x nvar x nshocks x nregimes [x npctg])
%==========================================================================
% [INF, SUP, MED, BAR] = TVARirband(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel.m
%   - VARopt : options (nsteps, ndraws, pctg, nreps_girf, etc.)
% -------------------------------------------------------------------------
% OUTPUT
%   - INF : lower band
%   - SUP : upper band
%   - MED : median
%   - BAR : mean
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------


%% Check inputs
%==========================================================================
if ~exist('TVAR','var')
    error('You need to provide TVAR structure, result of TVARmodel');
end

nsteps   = VARopt.nsteps;
ndraws   = VARopt.ndraws;
pctg     = VARopt.pctg;
nvar     = TVAR.nvar;
nlag     = TVAR.nlag;
nregimes = TVAR.nregimes;

% GIRF options
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

delay        = TVAR.delay;
thrvar_idx   = VARopt.thrvar_idx;
tartransform = VARopt.tartransform;


%% Dispatch by estimation method
%==========================================================================
if strcmp(TVAR.method, 'bayes')
    [INF, SUP, MED, BAR] = bands_bayes(TVAR, VARopt, nsteps, ndraws, pctg, ...
        nvar, nlag, nregimes, nreps_sim, shock_scale, delay, thrvar_idx, tartransform);
else
    [INF, SUP, MED, BAR] = bands_freq(TVAR, VARopt, nsteps, ndraws, pctg, ...
        nvar, nlag, nregimes, nreps_sim, shock_scale, delay, thrvar_idx, tartransform);
end

end


%% ========================================================================
%  BAYESIAN BANDS: GIRF at each posterior draw
%  ========================================================================
function [INF, SUP, MED, BAR] = bands_bayes(TVAR, VARopt, nsteps, ndraws, pctg, ...
    nvar, nlag, nregimes, nreps_sim, shock_scale, delay, thrvar_idx, tartransform)

fsize = TVAR.ndraws;
shock_pos = 1:nvar;

% Limit draws
nuse = min(fsize, ndraws);

% Pre-allocate: (nsteps x nvar x nshocks x nregimes x ndraws)
IR_all = nan(nsteps, nvar, nvar, nregimes, nuse);

disp('Computing Bayesian credible bands for TVAR GIRFs...');
for jj = 1:nuse
    if mod(jj, VARopt.mult) == 0
        fprintf('  Draw %d / %d\n', jj, nuse);
    end

    % Extract posterior draw (already in Mumtaz ordering)
    beta1 = TVAR.beta1_draws(jj,:)';
    beta2 = TVAR.beta2_draws(jj,:)';
    sigma1 = squeeze(TVAR.sigma1_draws(jj,:,:));
    sigma2 = squeeze(TVAR.sigma2_draws(jj,:,:));
    thresh = TVAR.thresh_draws(jj);

    % Check stability
    S1 = stability_tvar(beta1, nvar, nlag);
    S2 = stability_tvar(beta2, nvar, nlag);
    if S1 || S2
        continue;
    end

    % Check positive definiteness
    try
        chol(sigma1);
        chol(sigma2);
    catch
        continue;
    end

    % Compute GIRF at this posterior draw
    try
        GIRF = TVARgirf(beta1, sigma1, beta2, sigma2, thresh, ...
            TVAR.Y, TVAR.Ystar, nvar, nlag, delay, nsteps, nreps_sim, ...
            shock_pos, shock_scale, tartransform, thrvar_idx);

        IR_all(:,:,:,1,jj) = GIRF.regime1;
        IR_all(:,:,:,2,jj) = GIRF.regime2;
    catch
        continue;
    end
end
disp('-- Done!');

% Compute bands
[INF, SUP, MED, BAR] = compute_bands(IR_all, pctg, nsteps, nvar, nregimes);

end


%% ========================================================================
%  FREQUENTIST BANDS: bootstrap + GIRF at each draw
%  ========================================================================
function [INF, SUP, MED, BAR] = bands_freq(TVAR, VARopt, nsteps, ndraws, pctg, ...
    nvar, nlag, nregimes, nreps_sim, shock_scale, delay, thrvar_idx, tartransform)

const     = TVAR.const;
method    = VARopt.method;
ntotcoeff = TVAR.regime{1}.ntotcoeff;
shock_pos = 1:nvar;

% Stability threshold
maxEig_pt = max(cellfun(@(r) r.maxEig, TVAR.regime));
stab_thresh = max(1.05, maxEig_pt + 0.10);

% Pre-allocate
IR_all = nan(nsteps, nvar, nvar, nregimes, ndraws);

disp('Computing frequentist bootstrap bands for TVAR GIRFs...');
tt = 1;
ww = 1;
max_attempts = ndraws * 5;
attempts = 0;

while tt <= ndraws && attempts < max_attempts
    attempts = attempts + 1;

    if tt == VARopt.mult * ww
        disp(['  Loop ' num2str(tt) ' / ' num2str(ndraws) ' draws']);
        ww = ww + 1;
    end

    all_stable = true;
    beta1_boot = [];
    sigma1_boot = [];
    beta2_boot = [];
    sigma2_boot = [];

    for k = 1:nregimes
        reg = TVAR.regime{k};
        nobsk = reg.nobs;
        residk = reg.resid;

        % Resample residuals
        if strcmp(method, 'wild')
            rr = 1 - 2*(rand(nobsk,1) > 0.5);
            u = residk .* (rr * ones(1, nvar));
        else
            u = residk(ceil(nobsk * rand(nobsk, 1)), :);
        end

        % Generate artificial Y
        Y_art = reg.X * reg.Ft + u;

        % Re-estimate OLS
        Ft_draw = (reg.X' * reg.X) \ (reg.X' * Y_art);
        F_draw  = Ft_draw';
        resid_draw = Y_art - reg.X * Ft_draw;
        sigma_draw = (resid_draw' * resid_draw) / (nobsk - ntotcoeff);

        % Companion stability check
        Fcomp_draw = [F_draw(:, 1+const:nvar*nlag+const); ...
                      eye(nvar*(nlag-1)), zeros(nvar*(nlag-1), nvar)];
        maxEig_draw = max(abs(eig(Fcomp_draw)));

        if maxEig_draw >= stab_thresh
            all_stable = false;
            break;
        end

        % Convert to Mumtaz ordering [lags|const]
        beta_mumtaz = [Ft_draw(2:nvar*nlag+1, :); Ft_draw(1, :)];
        if k == 1
            beta1_boot = beta_mumtaz(:);
            sigma1_boot = sigma_draw;
        else
            beta2_boot = beta_mumtaz(:);
            sigma2_boot = sigma_draw;
        end
    end

    if ~all_stable
        continue;
    end

    % Check positive definiteness
    try
        chol(sigma1_boot);
        chol(sigma2_boot);
    catch
        continue;
    end

    % Compute GIRF at bootstrap draw
    try
        GIRF = TVARgirf(beta1_boot, sigma1_boot, beta2_boot, sigma2_boot, ...
            TVAR.thresh, TVAR.Y, TVAR.Ystar, nvar, nlag, delay, nsteps, ...
            nreps_sim, shock_pos, shock_scale, tartransform, thrvar_idx);

        IR_all(:,:,:,1,tt) = GIRF.regime1;
        IR_all(:,:,:,2,tt) = GIRF.regime2;
        tt = tt + 1;
    catch
        continue;
    end
end

if tt - 1 < ndraws
    warning('TVARirband: only %d of %d bootstrap draws were successful', tt-1, ndraws);
    IR_all = IR_all(:,:,:,:,1:tt-1);
end
disp('-- Done!');

% Compute bands
[INF, SUP, MED, BAR] = compute_bands(IR_all, pctg, nsteps, nvar, nregimes);

end


%% ========================================================================
%  HELPER: compute percentile bands from IR_all
%  ========================================================================
function [INF, SUP, MED, BAR] = compute_bands(IR_all, pctg, nsteps, nvar, nregimes)

npctg = numel(pctg);
draw_dim = 5;

% Drop NaN/complex draws
ndraws_actual = size(IR_all, draw_dim);
valid = true(1, ndraws_actual);
for jj = 1:ndraws_actual
    slice = IR_all(:,:,:,:,jj);
    if any(~isreal(slice(:))) || all(isnan(slice(:)))
        valid(jj) = false;
    end
end
IR_all = IR_all(:,:,:,:,valid);

if npctg > 1
    INF = zeros(nsteps, nvar, nvar, nregimes, npctg);
    SUP = zeros(nsteps, nvar, nvar, nregimes, npctg);
    for pp = 1:npctg
        INF(:,:,:,:,pp) = prctile(IR_all, (100-pctg(pp))/2, draw_dim);
        SUP(:,:,:,:,pp) = prctile(IR_all, 100-(100-pctg(pp))/2, draw_dim);
    end
else
    INF = prctile(IR_all, (100-pctg)/2, draw_dim);
    SUP = prctile(IR_all, 100-(100-pctg)/2, draw_dim);
end
MED = prctile(IR_all, 50, draw_dim);
BAR = nanmean(IR_all, draw_dim);

end
