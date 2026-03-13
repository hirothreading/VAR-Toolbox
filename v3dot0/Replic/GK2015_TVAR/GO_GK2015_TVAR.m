% Threshold VAR extension of Gertler and Karadi (2015, AEJ:M).
% Uses EBP as the threshold variable to identify regime-dependent
% monetary policy transmission.
%
% GIRFs (Generalized Impulse Response Functions) are computed following
% Koop, Pesaran & Potter (1996), allowing endogenous regime switching
% during the IRF horizon.
%==========================================================================
% The VAR Toolbox 3.0 is required to run this code. To get the
% latest version of the toolboxes visit:
% https://github.com/ambropo/VAR-Toolbox
%==========================================================================
% Ambrogio Cesa Bianchi
% ambrogio.cesabianchi@gmail.com


%% PRELIMINARIES
%==========================================================================
clear all; clear session; close all; clc
warning off all


%% LOAD DATA
%==========================================================================
% Load data (same dataset as the standard GK2015 replication)
[xlsdata, xlstext] = xlsread('../GK2015/GK2015_Data.xlsx','Sheet1');
dates = xlstext(3:end,1);
vnames_long = xlstext(1,2:end);
vnames = xlstext(2,2:end);
nvar = length(vnames);
data   = Num2NaN(xlsdata);
% Store variables in the structure DATA
for ii=1:length(vnames)
    DATA.(vnames{ii}) = data(:,ii);
end
% Convert the first date to numeric
year = str2double(xlstext{3,1}(1:4));
quarter = str2double(xlstext{3,1}(6));
% Observations
nobs = size(data,1);


%% SET UP TVAR
%==========================================================================
% Set endogenous variables (same as GK2015)
VARvnames_long = {'1yr rate';'CPI';'IP';'EBP'};
VARvnames      = {'gs1';'logcpi';'logip';'ebp'};
VARnvar        = length(VARvnames);

% VAR specification
VARnlags = 12; VARconst = 1;

% Create matrices of variables for the VAR
ENDO = nan(nobs,VARnvar);
for ii=1:VARnvar
    ENDO(:,ii) = DATA.(VARvnames{ii});
end


%% TVAR OPTIONS
%==========================================================================
TVARopt = TVARoption;
TVARopt.tvar_method = 'bayes';    % Bayesian Gibbs sampler
TVARopt.thrvar_idx  = 4;          % EBP as threshold variable
TVARopt.delay       = 1;          % 1-period lag of threshold variable
TVARopt.tartransform = 0;         % threshold in levels
TVARopt.nreps       = 20000;      % Gibbs replications
TVARopt.nburn       = 5000;       % burn-in
TVARopt.nskip       = 5;          % thinning
TVARopt.lambdaP     = 0.2;        % Minnesota prior: tightness on first lag
TVARopt.tauP        = 2.0;        % Minnesota prior: sum-of-coefficients
TVARopt.epsilonP    = 0.001;      % Minnesota prior: constant
TVARopt.tarvariance = 10;         % prior variance on threshold
TVARopt.tarscale    = 0.001;      % MH proposal variance
TVARopt.mult        = 1000;       % display frequency


%% ESTIMATE TVAR
%==========================================================================
disp('Estimating Threshold VAR (Bayesian)...')
[TVAR, TVARopt] = TVARmodel(ENDO, VARnlags, VARconst, TVARopt);

% Display estimation results
TVARopt.vnames = VARvnames_long;
TVARprint(TVAR, TVARopt);


%% CHOLESKY IDENTIFICATION: TVAR GIRFs
%==========================================================================
disp('Computing Cholesky GIRFs per regime...')
TVARopt.ident       = 'short';
TVARopt.nsteps      = 48;
TVARopt.nreps_girf  = 50;         % simulations per starting point
TVARopt.shock_scale = 1;          % 1 SD shock

% Point estimates
[IR_chol, TVAR] = TVARir(TVAR, TVARopt);

% Credible bands from posterior draws
TVARopt.ndraws = 200;
TVARopt.pctg   = 68;
[INF_chol, SUP_chol, MED_chol, BAR_chol] = TVARirband(TVAR, TVARopt);

% Plot: overlay mode (both regimes on same figure)
TVARopt.figname = 'GK2015_TVAR_';
TVARopt.FigSize = [20 16];
TVARopt.snames  = VARvnames_long;
TVARopt.regime_names = {'Low EBP','High EBP'};
TVARopt.tvar_plot_mode = 'overlay';
TVARirplot(IR_chol, TVARopt, INF_chol, SUP_chol);
disp('Cholesky GIRF plots saved.')


%% VARIANCE DECOMPOSITION (Cholesky only)
%==========================================================================
disp('Computing variance decomposition per regime...')
TVARopt.ident = 'short';
TVARopt.figname = 'GK2015_TVAR_';
[VD, TVAR] = TVARvd(TVAR, TVARopt);
TVARvdplot(VD, TVARopt);
disp('VD plots saved.')


%% HISTORICAL DECOMPOSITION (Cholesky only)
%==========================================================================
disp('Computing historical decomposition...')
TVARopt.ident = 'short';
TVARopt.figname = 'GK2015_TVAR_';
[HD, TVAR] = TVARhd(TVAR, TVARopt);
TVARhdplot(HD, TVARopt, TVAR);
disp('HD plots saved.')


disp(' ')
disp('=== GK2015 TVAR Replication Complete ===')
