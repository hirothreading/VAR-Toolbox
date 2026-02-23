function TVARopt = TVARoption
%==========================================================================
% Default options for Threshold VAR (TVAR) analysis. Inherits all fields
% from VARoption and appends TVAR-specific settings.
%==========================================================================
% TVARopt = TVARoption
%==========================================================================
% VAR Toolbox 3.0 - TVAR Extension
%==========================================================================

% Inherit all standard VAR options
TVARopt = VARoption;

% -------------------------------------------------------------------
% Threshold variable
% -------------------------------------------------------------------
TVARopt.thrvar_idx = 1;       % Column index of threshold variable in ENDO.
                               % Set to 0 when using an external threshold
                               % variable (passed as THRVAR_EX to TVARmodel).

TVARopt.delay      = 1;       % Delay d: threshold variable at time t is
                               % q_{t-d}. Must satisfy 1 <= delay <= nlag.

% -------------------------------------------------------------------
% Threshold search
% -------------------------------------------------------------------
TVARopt.nthresh    = 1;       % Number of thresholds (only 1 supported).
                               % Two regimes: q_{t-d} <= thresh (regime 1)
                               % and q_{t-d} > thresh (regime 2).

TVARopt.trim       = 0.15;    % Fraction of obs trimmed from each tail of
                               % the threshold variable when defining the
                               % grid. Ensures each regime has enough obs.
                               % Default 0.15 follows Hansen (1999).

TVARopt.ngrid      = 300;     % Number of equally-spaced grid points
                               % searched for the threshold value.

% -------------------------------------------------------------------
% Display / plotting
% -------------------------------------------------------------------
TVARopt.rnames     = {'Regime 1', 'Regime 2'};
                               % Regime labels used in TVARprint and
                               % TVARirplot. Override with meaningful names,
                               % e.g. {'Recession', 'Expansion'}.
