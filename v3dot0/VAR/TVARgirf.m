function [GIRF] = TVARgirf(beta1, sigma1, beta2, sigma2, thresh, ...
    Y, Ystar, nvar, nlag, delay, nsteps, nreps_sim, shock_pos, shock_scale, ...
    tartransform, thrvar_idx)
%==========================================================================
% Compute Generalized Impulse Response Functions (GIRFs) for a two-regime
% Threshold VAR, following the Koop, Pesaran & Potter (1996) approach as
% implemented in Haroon Mumtaz's TVAR Toolkit.
%
% The GIRF is defined as:
%   GIRF(h) = E[y_{t+h} | shock, history] - E[y_{t+h} | history]
%
% For each historical starting point in a given regime, the function:
%   1. Simulates an unshocked path forward (with random N(0,1) innovations)
%   2. Simulates a shocked path (same innovations + structural shock at h=0)
%   3. Takes the difference as the IRF for that starting point
%   4. Averages across simulations and starting points
%
% The system can switch regimes endogenously during the IRF horizon.
%==========================================================================
% [GIRF] = TVARgirf(beta1, sigma1, beta2, sigma2, thresh, ...
%     Y, Ystar, nvar, nlag, delay, nsteps, nreps_sim, shock_pos, ...
%     shock_scale, tartransform, thrvar_idx)
% -------------------------------------------------------------------------
% INPUT
%   - beta1      : regime 1 coefficients, (nvar*(nvar*nlag+1) x 1) in
%                  Mumtaz ordering [lag1 lag2 ... lagL constant]
%   - sigma1     : regime 1 covariance matrix (nvar x nvar)
%   - beta2      : regime 2 coefficients (same format)
%   - sigma2     : regime 2 covariance matrix (nvar x nvar)
%   - thresh     : threshold value (scalar)
%   - Y          : data matrix used for estimation (nobse x nvar)
%   - Ystar      : threshold variable (nobse x 1), already transformed/lagged
%   - nvar       : number of endogenous variables
%   - nlag       : number of lags
%   - delay      : delay for threshold variable
%   - nsteps     : IRF horizon
%   - nreps_sim  : number of simulation repetitions per starting point
%   - shock_pos  : position (variable index) of the shock (scalar or vector)
%   - shock_scale: size of shock in standard deviations (scalar)
%   - tartransform: 0 = level, 1 = annual growth (4-period diff)
%   - thrvar_idx : column index of threshold variable in Y
% -------------------------------------------------------------------------
% OUTPUT
%   - GIRF : structure with fields:
%       .regime1 : (nsteps x nvar x nshocks) GIRFs conditioned on regime 1
%       .regime2 : (nsteps x nvar x nshocks) GIRFs conditioned on regime 2
% =========================================================================
% VAR Toolbox 3.0 — TVAR Extension
% Adapted from Haroon Mumtaz's getimpulse.m
% -------------------------------------------------------------------------


%% Setup
%==========================================================================
ncoeff = nvar * nlag + 1;  % Mumtaz ordering: [lags | constant]
b1x = reshape(beta1, ncoeff, nvar);
b2x = reshape(beta2, ncoeff, nvar);
csigma1 = chol(sigma1);  % upper triangular Cholesky
csigma2 = chol(sigma2);

% Transformation lag
if tartransform == 1
    Lx = 4;
else
    Lx = 0;
end
LL = max([nlag, delay, Lx + delay]);

% Determine number of shocks
if isscalar(shock_pos)
    nshocks = 1;
    shock_positions = shock_pos;
else
    nshocks = length(shock_pos);
    shock_positions = shock_pos;
end

% Regime classification
e1 = (Ystar <= thresh);
T = size(Y, 1);

% Apply transformation to Y for threshold evaluation during simulation
if tartransform == 1
    Ym = [nan(Lx, nvar); Y(Lx+1:end,:) - Y(1:end-Lx,:)];
else
    Ym = Y;
end

%% Build histories for each regime
%==========================================================================
% History = the LL most recent observations ending at time t
history = cell(T, 1);
historym = cell(T, 1);
for j = LL:T
    history{j} = Y(j-LL+1:j, :);
    historym{j} = Ym(j-LL+1:j, :);
end

% Separate histories by regime
T1 = sum(e1(LL:end));  % regime 1 observations with sufficient history
T2 = sum(~e1(LL:end));

history1 = cell(T1, 1);  historym1 = cell(T1, 1);
history2 = cell(T2, 1);  historym2 = cell(T2, 1);
jj1 = 1;  jj2 = 1;
for j = LL:T
    if e1(j) == 1
        history1{jj1} = history{j};
        historym1{jj1} = historym{j};
        jj1 = jj1 + 1;
    else
        history2{jj2} = history{j};
        historym2{jj2} = historym{j};
        jj2 = jj2 + 1;
    end
end


%% Compute GIRFs for each regime
%==========================================================================
GIRF.regime1 = compute_girf_regime(history1, historym1, T1, LL, ...
    nvar, nlag, nsteps, nreps_sim, b1x, b2x, csigma1, csigma2, ...
    sigma1, sigma2, thresh, thrvar_idx, delay, Lx, tartransform, ...
    shock_positions, shock_scale, nshocks);

GIRF.regime2 = compute_girf_regime(history2, historym2, T2, LL, ...
    nvar, nlag, nsteps, nreps_sim, b1x, b2x, csigma1, csigma2, ...
    sigma1, sigma2, thresh, thrvar_idx, delay, Lx, tartransform, ...
    shock_positions, shock_scale, nshocks);

end


%% ========================================================================
%  Helper: compute GIRF for one regime's starting points
%  ========================================================================
function girf_out = compute_girf_regime(histories, historiesm, Tk, LL, ...
    nvar, nlag, nsteps, nreps_sim, b1x, b2x, csigma1, csigma2, ...
    sigma1, sigma2, thresh, thrvar_idx, delay, Lx, tartransform, ...
    shock_positions, shock_scale, nshocks)

girf_out = zeros(nsteps, nvar, nshocks);

% Count valid starting points
nvalid = 0;
for t = 1:Tk
    if ~isempty(histories{t})
        nvalid = nvalid + 1;
    end
end

if nvalid == 0
    warning('TVARgirf: no valid starting points for this regime');
    return;
end

% Cholesky factors for structural identification
A01 = chol(sigma1);  % upper triangular
A02 = chol(sigma2);

% Loop over starting points
girf_accum = zeros(nsteps, nvar, nshocks);

for t = 1:Tk
    Y0 = histories{t};
    Y0m = historiesm{t};
    if isempty(Y0); continue; end

    % Loop over shocks
    for ss = 1:nshocks
        pos = shock_positions(ss);
        ir_t = simulate_girf(Y0, Y0m, nvar, nlag, nsteps, nreps_sim, ...
            b1x, b2x, csigma1, csigma2, A01, A02, thresh, ...
            thrvar_idx, delay, Lx, tartransform, pos, shock_scale, LL);
        girf_accum(:,:,ss) = girf_accum(:,:,ss) + ir_t;
    end
end

girf_out = girf_accum / nvalid;

end


%% ========================================================================
%  Core simulation: GIRF for one starting point and one shock
%  ========================================================================
function ir = simulate_girf(Y0, Y0m, nvar, nlag, nsteps, nreps_sim, ...
    b1x, b2x, csigma1, csigma2, A01, A02, thresh, ...
    thrvar_idx, delay, Lx, tartransform, pos, shock_scale, LL)

% Accumulate shocked and unshocked paths
yy  = zeros(nsteps + LL, nvar);   % unshocked path accumulator
yy1 = zeros(nsteps + LL, nvar);   % shocked path accumulator

for ii = 1:nreps_sim
    % Initialize paths with history
    yhat  = zeros(nsteps + LL, nvar);
    yhat1 = zeros(nsteps + LL, nvar);
    yhat(1:LL, :)  = Y0;
    yhat1(1:LL, :) = Y0;

    % Transformed paths for threshold evaluation
    yhatm  = zeros(nsteps + LL, nvar);
    yhat1m = zeros(nsteps + LL, nvar);
    yhatm(1:LL, :)  = Y0m;
    yhat1m(1:LL, :) = Y0m;

    % Threshold variable paths
    ystar  = zeros(nsteps + LL, 1);
    ystar1 = zeros(nsteps + LL, 1);

    % Simulate forward
    for fi = LL+1 : nsteps+LL
        % Build regressor vectors [lag1 lag2 ... lagL constant]
        xhat  = ones(1, nvar*nlag + 1);
        xhat1 = ones(1, nvar*nlag + 1);
        jix = 1;
        for ji = 1:nlag
            xhat(jix:jix+nvar-1)  = yhat(fi-ji, :);
            xhat1(jix:jix+nvar-1) = yhat1(fi-ji, :);
            jix = jix + nvar;
        end

        % Update transformed values for threshold evaluation
        if tartransform == 0
            yhatm(fi,:)  = yhat(fi-1,:);   % will be overwritten below
            yhat1m(fi,:) = yhat1(fi-1,:);
        elseif tartransform == 1
            if fi-delay >= 1 && fi-Lx-delay >= 1
                yhatm(fi,:)  = yhat(fi-delay,:)  - yhat(fi-Lx-delay,:);
                yhat1m(fi,:) = yhat1(fi-delay,:) - yhat1(fi-Lx-delay,:);
            end
        end

        % Evaluate threshold variable
        if tartransform == 0
            ystar(fi)  = yhat(fi-delay, thrvar_idx);
            ystar1(fi) = yhat1(fi-delay, thrvar_idx);
        else
            ystar(fi)  = yhatm(fi, thrvar_idx);
            ystar1(fi) = yhat1m(fi, thrvar_idx);
        end

        % Determine regime for each path
        e1  = (ystar(fi)  <= thresh);
        e1s = (ystar1(fi) <= thresh);

        % Simulate
        if fi == LL + 1
            % First period: no random innovation for unshocked path,
            % structural shock for shocked path
            uu = zeros(1, nvar);
            yhat(fi,:)  = (xhat  * b1x + uu * csigma1) * e1 + ...
                          (xhat  * b2x + uu * csigma2) * (~e1);

            uu_shock = zeros(1, nvar);
            uu_shock(pos) = shock_scale;
            yhat1(fi,:) = (xhat1 * b1x + uu_shock * A01) * e1s + ...
                          (xhat1 * b2x + uu_shock * A02) * (~e1s);
        else
            % Subsequent periods: same random innovations for both paths
            uu = randn(1, nvar);
            yhat(fi,:)  = (xhat  * b1x + uu * csigma1) * e1 + ...
                          (xhat  * b2x + uu * csigma2) * (~e1);
            yhat1(fi,:) = (xhat1 * b1x + uu * csigma1) * e1s + ...
                          (xhat1 * b2x + uu * csigma2) * (~e1s);
        end

        % Update transformed value at current period (for tartransform==0)
        if tartransform == 0
            yhatm(fi,:)  = yhat(fi,:);
            yhat1m(fi,:) = yhat1(fi,:);
        end
    end

    yy  = yy  + yhat;
    yy1 = yy1 + yhat1;
end

% Average across simulations
yy  = yy  / nreps_sim;
yy1 = yy1 / nreps_sim;

% GIRF = difference between shocked and unshocked paths
ir = yy1(LL+1:end, :) - yy(LL+1:end, :);

end
