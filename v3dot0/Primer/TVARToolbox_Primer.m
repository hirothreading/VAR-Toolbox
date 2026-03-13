% TVAR Toolbox Primer
% =======================================================================
% This script demonstrates the Threshold VAR (TVAR) capabilities of the
% VAR Toolbox 3.0. It covers:
%   1. Data preparation
%   2. TVARtest — Hansen (1999) linearity test
%   3. TVARmodel — TVAR estimation (frequentist and Bayesian)
%   4. TVARprint — Display results
%   5. TVARir / TVARirband / TVARirplot — Generalized IRFs (KPP 1996)
%   6. TVARvd / TVARvdplot — Variance decomposition
%   7. TVARhd / TVARhdplot — Historical decomposition
% =======================================================================
% VAR Toolbox 3.0
% Ambrogio Cesa-Bianchi
% ambrogio.cesabianchi@gmail.com
% =======================================================================


%% PRELIMINARIES
%==========================================================================
clear all; clear session; close all; clc
warning off all

% Add the toolbox to the path (adjust as needed)
% addpath(genpath('/path/to/VAR-Toolbox/v3dot0/'))


%% 1. LOAD DATA
%==========================================================================
% Use the Simple_Data dataset for illustration
[xlsdata, xlstext] = xlsread('data/Simple_Data.xlsx','Sheet1');
dates = xlstext(3:end,1);
vnames_long = xlstext(1,2:end);
vnames = xlstext(2,2:end);
nvar = length(vnames);
data = Num2NaN(xlsdata);

% Store variables
for ii=1:length(vnames)
    DATA.(vnames{ii}) = data(:,ii);
end

nobs = size(data,1);
disp(['Loaded ' num2str(nobs) ' observations, ' num2str(nvar) ' variables.'])


%% 2. SET UP ENDOGENOUS VARIABLES
%==========================================================================
% Select variables for the TVAR
VARvnames_long = vnames_long(1:3);  % first 3 variables
VARvnames      = vnames(1:3);
VARnvar        = length(VARvnames);

% Build ENDO matrix
ENDO = nan(nobs, VARnvar);
for ii = 1:VARnvar
    ENDO(:,ii) = DATA.(VARvnames{ii});
end

% VAR specification
VARnlags = 4;
VARconst = 1;


%% 3. TVAR OPTIONS
%==========================================================================
TVARopt = TVARoption;

% Threshold variable: use the first endogenous variable
TVARopt.thrvar_idx   = 1;        % column index in ENDO
TVARopt.delay        = 1;        % 1-period lag
TVARopt.tartransform = 0;        % threshold in levels
TVARopt.tvar_method  = 'freq';   % frequentist estimation
TVARopt.trim         = 0.15;     % 15% trimming
TVARopt.ngrid        = 300;      % grid points


%% 4. LINEARITY TEST (Hansen 1999)
%==========================================================================
% Test H0: linear VAR vs H1: two-regime TVAR
% (This can be slow — reduce ndraws for a quick check)
disp('Running linearity test (this may take a few minutes)...')
TVARopt.ndraws = 100;    % use 500+ for publication
TVARopt.method = 'wild';
test = TVARtest(ENDO, VARnlags, VARconst, TVARopt);
disp(['Linearity test p-value: ' num2str(test.pval,'%.4f')])


%% 5. ESTIMATE TVAR
%==========================================================================
disp('Estimating TVAR...')
[TVAR, TVARopt] = TVARmodel(ENDO, VARnlags, VARconst, TVARopt);

% Add variable labels
TVARopt.vnames = VARvnames_long;
TVARopt.snames = VARvnames_long;


%% 6. DISPLAY RESULTS
%==========================================================================
TVARprint(TVAR, TVARopt);


%% 7. CHOLESKY GENERALIZED IMPULSE RESPONSES
%==========================================================================
% GIRFs are computed following Koop, Pesaran & Potter (1996).
% The system is simulated forward from historical starting points,
% allowing endogenous regime switching during the IRF horizon.
disp('Computing Cholesky GIRFs...')
TVARopt.ident       = 'short';
TVARopt.nsteps      = 24;
TVARopt.nreps_girf  = 50;     % simulations per starting point
TVARopt.shock_scale = 1;      % 1 SD shock
TVARopt.method      = 'wild';
TVARopt.ndraws      = 200;
TVARopt.pctg        = 68;
TVARopt.FigSize     = [18 14];
TVARopt.quality     = 1;
TVARopt.figname     = 'graphics/TVAR_chol_';
TVARopt.regime_names = {'Regime 1','Regime 2'};

% Point estimates
[IR, TVAR] = TVARir(TVAR, TVARopt);

% Confidence bands
[INF, SUP, MED, BAR] = TVARirband(TVAR, TVARopt);

% Plot — overlay mode (both regimes on same figure)
TVARopt.tvar_plot_mode = 'overlay';
TVARirplot(IR, TVARopt, INF, SUP);
disp('Cholesky GIRF plots saved to graphics/')


%% 8. VARIANCE DECOMPOSITION
%==========================================================================
disp('Computing variance decomposition...')
TVARopt.figname = 'graphics/TVAR_chol_';
[VD, TVAR] = TVARvd(TVAR, TVARopt);
TVARvdplot(VD, TVARopt);
disp('VD plots saved to graphics/')


%% 9. HISTORICAL DECOMPOSITION
%==========================================================================
disp('Computing historical decomposition...')
TVARopt.figname = 'graphics/TVAR_chol_';
[HD, TVAR] = TVARhd(TVAR, TVARopt);
TVARhdplot(HD, TVARopt, TVAR);
disp('HD plots saved to graphics/')


disp(' ')
disp('=== TVAR Toolbox Primer Complete ===')
