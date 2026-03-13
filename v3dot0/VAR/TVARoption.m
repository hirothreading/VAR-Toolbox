function TVARopt = TVARoption
%==========================================================================
% Default options for Threshold VAR (TVAR) analysis. This function is
% called automatically by TVARmodel. Users can override any field after
% calling it.
%==========================================================================
% TVARopt = TVARoption
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% Ambrogio Cesa-Bianchi
% -------------------------------------------------------------------------

%% Inherit all standard VAR options
% -------------------------------------------------------------------------
TVARopt = VARoption;

%% TVAR estimation method
% -------------------------------------------------------------------------
TVARopt.tvar_method  = 'freq';   % 'freq' (Hansen grid search + OLS) or 'bayes' (Gibbs sampler)

%% Threshold specification
% -------------------------------------------------------------------------
TVARopt.thrvar_idx   = 1;        % column index of threshold variable in ENDO (0 = external)
TVARopt.delay        = 1;        % delay d for threshold variable: q_{t-d}
TVARopt.tartransform = 0;        % transformation of threshold variable: 0=level, 1=annual growth (4-period diff)

%% Frequentist options (method='freq')
% -------------------------------------------------------------------------
TVARopt.trim         = 0.15;     % trimming fraction: exclude bottom/top quantiles from grid search
TVARopt.ngrid        = 300;      % number of grid points for threshold search

%% Bayesian options (method='bayes')
% -------------------------------------------------------------------------
TVARopt.nreps        = 20000;    % total Gibbs sampler replications
TVARopt.nburn        = 5000;     % burn-in period
TVARopt.nskip        = 5;        % thinning: save every nskip-th draw
TVARopt.lambdaP      = 0.2;      % Minnesota prior: tightness on first lag coefficients
TVARopt.tauP         = 2.0;      % Minnesota prior: tightness on sum-of-coefficients
TVARopt.epsilonP     = 0.001;    % Minnesota prior: tightness on constant
TVARopt.tarvariance  = 10;       % prior variance on threshold value (mean = sample mean of Ystar)
TVARopt.tarscale     = 0.001;    % random walk MH proposal variance for threshold
TVARopt.maxtrys      = 1000;     % max draws to find a stable coefficient vector

%% GIRF options (Generalized Impulse Response Functions)
% -------------------------------------------------------------------------
TVARopt.nreps_girf     = 50;       % number of simulation repetitions per starting point
TVARopt.shock_scale    = 1;        % shock size in standard deviations

%% Plotting options
% -------------------------------------------------------------------------
TVARopt.tvar_plot_mode = 'overlay'; % 'overlay' (both regimes on one figure) or 'separate' (one figure per regime)
TVARopt.regime_names   = {'Regime 1', 'Regime 2'}; % labels for the two regimes
