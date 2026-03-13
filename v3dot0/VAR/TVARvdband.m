function [INF, SUP, MED, BAR] = TVARvdband(TVAR, VARopt)
%==========================================================================
% Compute confidence / credible bands for TVAR variance decompositions.
% Uses the same bootstrap / posterior-draw approach as TVARirband.
%
% Not available for 'iv' identification.
%==========================================================================
% [INF, SUP, MED, BAR] = TVARvdband(TVAR, VARopt)
% -------------------------------------------------------------------------
% INPUT
%   - TVAR   : structure from TVARmodel.m
%   - VARopt : options (nsteps, ndraws, pctg, method, ident)
% -------------------------------------------------------------------------
% OUTPUT
%   - INF, SUP, MED, BAR : (nsteps x nshocks x nvar x nregimes [x npctg])
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% -------------------------------------------------------------------------

if strcmp(VARopt.ident, 'iv')
    error('TVARvdband: FEVD bands not available with iv identification');
end

nsteps   = VARopt.nsteps;
ndraws   = VARopt.ndraws;
pctg     = VARopt.pctg;
nvar     = TVAR.nvar;
nregimes = TVAR.nregimes;

if strcmp(TVAR.method, 'bayes')
    [INF, SUP, MED, BAR] = vdband_bayes(TVAR, VARopt, nsteps, ndraws, pctg, nvar, nregimes);
else
    [INF, SUP, MED, BAR] = vdband_freq(TVAR, VARopt, nsteps, ndraws, pctg, nvar, nregimes);
end
end


%% Bayesian bands
function [INF, SUP, MED, BAR] = vdband_bayes(TVAR, VARopt, nsteps, ndraws, pctg, nvar, nregimes)

fsize = TVAR.ndraws;
nlag  = TVAR.nlag;
const = TVAR.const;
ncoeff_bayes = nvar * nlag + 1;
nuse = min(fsize, ndraws);

VD_all = nan(nsteps, nvar, nvar, nregimes, nuse);

disp('Computing Bayesian VD bands...');
for jj = 1:nuse
    if mod(jj, VARopt.mult) == 0
        fprintf('  Draw %d / %d\n', jj, nuse);
    end
    for k = 1:nregimes
        if k == 1
            beta_draw = TVAR.beta1_draws(jj,:)';
            sigma_draw = squeeze(TVAR.sigma1_draws(jj,:,:));
        else
            beta_draw = TVAR.beta2_draws(jj,:)';
            sigma_draw = squeeze(TVAR.sigma2_draws(jj,:,:));
        end
        beta_mat = reshape(beta_draw, ncoeff_bayes, nvar);
        Ftk = [beta_mat(end,:); beta_mat(1:nvar*nlag,:)];
        Fk = Ftk';
        Fcomp_k = [Fk(:, 1+const:nvar*nlag+const); ...
                   eye(nvar*(nlag-1)), zeros(nvar*(nlag-1), nvar)];
        if max(abs(eig(Fcomp_k))) >= 0.9999; continue; end

        VAR_draw = TVAR.regime{k};
        VAR_draw.Ft = Ftk; VAR_draw.F = Fk;
        VAR_draw.sigma = sigma_draw; VAR_draw.Fcomp = Fcomp_k;
        VAR_draw.B = [];
        try
            [VD_draw, ~] = VARvd(VAR_draw, VARopt);
            VD_all(:,:,:,k,jj) = VD_draw;
        catch; end
    end
end
disp('-- Done!');
[INF, SUP, MED, BAR] = compute_vd_bands(VD_all, pctg);
end


%% Frequentist bands (fixed-regime bootstrap)
function [INF, SUP, MED, BAR] = vdband_freq(TVAR, VARopt, nsteps, ndraws, pctg, nvar, nregimes)

method = VARopt.method;
const  = TVAR.const;
ntotcoeff = TVAR.regime{1}.ntotcoeff;
nlag   = TVAR.nlag;

% Stability threshold: at least as large as point estimate's max eigenvalue
maxEig_pt = max(cellfun(@(r) r.maxEig, TVAR.regime));
stab_thresh = max(1.05, maxEig_pt + 0.10);

VD_all = nan(nsteps, nvar, nvar, nregimes, ndraws);

disp('Computing frequentist VD bands...');
tt = 1; ww = 1;
while tt <= ndraws
    if tt == VARopt.mult * ww
        disp(['  Loop ' num2str(tt) ' / ' num2str(ndraws)]); ww = ww+1;
    end
    TVAR_draw = TVAR;
    ok = true;
    for k = 1:nregimes
        reg = TVAR.regime{k};
        nobsk = reg.nobs;
        if strcmp(method,'bs')
            u = reg.resid(ceil(nobsk*rand(nobsk,1)),:);
        else
            rr = 1-2*(rand(nobsk,1)>0.5);
            u = reg.resid .* (rr*ones(1,nvar));
        end
        Y_art = reg.X * reg.Ft + u;
        Ft_d = (reg.X'*reg.X)\(reg.X'*Y_art);
        F_d = Ft_d';
        resid_d = Y_art - reg.X*Ft_d;
        sigma_d = (resid_d'*resid_d)/(nobsk - ntotcoeff);
        Fcomp_d = [F_d(:,1+const:nvar*nlag+const); eye(nvar*(nlag-1)), zeros(nvar*(nlag-1),nvar)];
        if max(abs(eig(Fcomp_d))) >= stab_thresh; ok = false; break; end
        TVAR_draw.regime{k}.Ft = Ft_d; TVAR_draw.regime{k}.F = F_d;
        TVAR_draw.regime{k}.sigma = sigma_d; TVAR_draw.regime{k}.Fcomp = Fcomp_d;
        TVAR_draw.regime{k}.B = [];
    end
    if ~ok; continue; end
    try
        [VD_draw, ~] = TVARvd(TVAR_draw, VARopt);
        VD_all(:,:,:,:,tt) = VD_draw; tt = tt+1;
    catch; end
end
disp('-- Done!');
[INF, SUP, MED, BAR] = compute_vd_bands(VD_all, pctg);
end


%% Helper
function [INF, SUP, MED, BAR] = compute_vd_bands(VD_all, pctg)
npctg = numel(pctg);
dd = 5;
if npctg > 1
    nsteps = size(VD_all,1); nvar = size(VD_all,2); nregimes = size(VD_all,4);
    INF = zeros(nsteps,nvar,nvar,nregimes,npctg);
    SUP = zeros(nsteps,nvar,nvar,nregimes,npctg);
    for pp = 1:npctg
        INF(:,:,:,:,pp) = prctile(VD_all,(100-pctg(pp))/2,dd);
        SUP(:,:,:,:,pp) = prctile(VD_all,100-(100-pctg(pp))/2,dd);
    end
else
    INF = prctile(VD_all,(100-pctg)/2,dd);
    SUP = prctile(VD_all,100-(100-pctg)/2,dd);
end
MED = prctile(VD_all,50,dd);
BAR = nanmean(VD_all,dd);
end
