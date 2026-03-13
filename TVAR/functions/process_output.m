% process_output:
% The script loads the output of a simulation and gets log-scores, PITs and
% other stuff. This version: P Alessandri, October 2012
% -------------------------------------------------------------------------

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET-UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set (global) parameters for functions below.
% NOTE: these must match those used in the simulation to be analysed.

global strFolder minfile maxfile REPS BURN HORZ T0 T00 VarBench VarAsset VarBank
 
% Bank of Italy (Osiride) setting:
% strFolder = '/home/user/wire141/private/september2012' ;
% strFolder = 'D:\My_Documents\Density_forecasts\Code\september2012';
strFolder = 'D:\My_Documents\Density_forecasts\Code\TVAR\files\';
strFolder = 'd:\dati\Density forecasts\Code\TVAR\files';

% Folder with forecasts (could be just 'forecasts'):
strFolderFc = 'forecasts';

minfile=1; % 1st file used to be indexed 0 [PG 7.1.13]
maxfile=355; 

REPS    = 20000;
BURN    = 15000;
HORZ    = 12;
T0      = [];   % Not used [PG 7.1.13]
T00     = 120;

% Choose 'target' variable:
iii = 1;
% This is the variable for which predictive densities will be analysed. The
% variable must appear in the same position in every model.

strModels = {'benchmark', 'TAR', '-'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get data:
matDataBench = get_data('data', 'benchmark');
vecGdp = matDataBench(:, 1);

% Get simulated values for all models (various horizons, h):
matGdpM1_h1 = extract_paths(strFolderFc, 'benchmark', iii, 1);
matGdpM2_h1  = extract_paths(strFolderFc, 'TAR', iii, 1);

matGdpM1_h6 = extract_paths(strFolderFc, 'benchmark', iii, 6);
matGdpM2_h6 = extract_paths(strFolderFc, 'TAR', iii, 6);

matGdpM1_h12 = extract_paths(strFolderFc, 'benchmark', iii, 12);
matGdpM2_h12 = extract_paths(strFolderFc, 'TAR', iii, 12);
% NOTE: PITs and Giacomini-White test must be based on 1-step-ahead
% forecasts.

% Get simulated 1Y cumulative outturn for all models:
% % matGdpBench_cum4 = generate_cumdata(matGdpBench_1, 4) / 4;
% % matGdpAsset_cum4 = generate_cumdata(matGdpAsset_1, 4) / 4;
% % matGdpBank_cum4  = generate_cumdata(matGdpBank_1, 4)  / 4;
% % CRAP! Cannot get 1Y forecast by adding up 4 1Q forecasts done at
% different points in time! Instead, rewrite extract_paths with a "type"
% input: type=normal for current use, type=cumulative to get cum forecast

% Get cumulative 1Y observed output growth:
vecGdp_cum6  = generate_cumdata(vecGdp', 6)' / 6;
vecGdp_cum12 = generate_cumdata(vecGdp', 12)' / 12;
% Note we divide by the number of periods because GDP growth is at annual 
% rates in the data. Transpositions needed to keep data vector as a column.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORECAST EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RMSE = struct;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DENSITY EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Log scores (LS)
% -------------------
LS = struct;
LS.M1 = [  get_logscore(matGdpM1_h1, vecGdp, 1);
           get_logscore(matGdpM1_h6, vecGdp, 6);
           get_logscore(matGdpM1_h12, vecGdp, 12)];

LS.M2 = [  get_logscore(matGdpM2_h1, vecGdp, 1);
           get_logscore(matGdpM2_h6, vecGdp, 6);
           get_logscore(matGdpM2_h12, vecGdp, 12)];

LS.M3 = NaN * LS.M1;

LS.averages = [nanmean(LS.M1, 2)  nanmean(LS.M2, 2) nanmean(LS.M3, 2)];
% Each LS.(model) matrix is k*T, where k is the number of horizons
% considered and T is the number of forecasts. LS.averages is k*m matrix
% with the average log-score for each model (m) at each horizon (k).

% Weighted log-scores:
% ---------------------
% Get the Amisano-Giacomini (2007) weights:
vecWeightsLT = get_LS_weights(vecGdp, 'recursive', 'ltail');
vecWeightsTs = get_LS_weights(vecGdp, 'recursive', 'tails');

% Get weighted log-scores (LEFT TAIL)
matWlsLt = NaN(3,3);
dblM1h1  = nansum( LS.M1(1,:) .* vecWeightsLT' ) / sum(~isnan(LS.M1(1,:)));
dblM2h1  = nansum( LS.M2(1,:) .* vecWeightsLT' ) / sum(~isnan(LS.M2(1,:)));
dblM3h1  = NaN;
matWlsLt(1, :) = [dblM1h1 dblM2h1 dblM3h1];

% Get weighted log-scores (BOTH TAILS)
matWlsTs = NaN(3,3);
dblM1h1  = nansum( LS.M1(1,:) .* vecWeightsTs' ) / sum(~isnan(LS.M1(1,:)));
dblM2h1  = nansum( LS.M2(1,:) .* vecWeightsTs' ) / sum(~isnan(LS.M2(1,:)));
dblM3h1  = NaN;
matWlsTs(1, :) = [dblM1h1 dblM2h1 dblM3h1];
% ----------------------
% NOTE: we (1) use nansum to ignore missing LS values; and (2) divide by the
% number of actual available log-scores to get the mean (using T0 or
% maxfile is wrong, because some forecasts cannot be scored).


% 2. Prob Integral Transforms (PIT)
PIT    = struct;
PIT.M1 = get_PITs(matGdpM1_h1, vecGdp); 
PIT.M2 = get_PITs(matGdpM2_h1, vecGdp); 
PIT.M3 = []; 

% 3. Inverse-normal transform of PITs:
PITin    = struct;
PITin.M1 = icdf('Normal', PIT.M1, 0, 1);
PITin.M2 = icdf('Normal', PIT.M2, 0, 1);
PITin.M3 = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHARTS etc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shoot average log-scores on screen:
disp('Average log-scores')
disp('------------------------------------')
disp(strModels)
disp(LS.averages)

% Shoot average weighted log-scores on screen:
disp('Average weighted log-scores, left tail')
disp('------------------------------------')
disp(strModels)
disp(matWlsLt)

disp('Average weighted log-scores, both tails')
disp('------------------------------------')
disp(strModels)
disp(matWlsTs)


% Histograms of PIT:
figure('Name', 'PITs')
kk = length(PIT.M1)/20;
subplot(1, 3, 1), hist(PIT.M1, kk), axis square, title(strModels(1))
subplot(1, 3, 2), hist(PIT.M2, kk), axis square, title(strModels(2))
subplot(1, 3, 3), hist(PIT.M3, kk), axis square, title(strModels(3))

% Histograms of inv-normal PITs:
figure('Name', 'PITs - inverse N')
subplot(1, 3, 1), hist(PITin.M1, kk), axis square, title(strModels(1))
subplot(1, 3, 2), hist(PITin.M2, kk), axis square, title(strModels(2))
subplot(1, 3, 3), hist(PITin.M3, kk), axis square, title(strModels(3))


% PLot log-scores over time 
figure('Name', 'Log-scores, 1Q ahead horizon')
title('Log-scores for 1Q-ahead predictions')
hold on
plot(LS.M1(1,:), 'b')
plot(LS.M2(1,:), 'r')
plot(LS.M3(1,:), 'k--')
hold off
xlabel('Time')
legend(strModels)

% PLot log-scores over time, various horizons
figure('Name', 'Log-scores')
title('Log scores')
subplot(2, 2, 1)
    hold on
    title('1 step ahead')
    plot(LS.M1(1,:), 'b')
    plot(LS.M2(1,:), 'r')
    plot(LS.M3(1,:), 'k--')
subplot(2, 2, 2)
    hold on
    title('4 steps ahead')
    plot(LS.M1(2,:), 'b')
    plot(LS.M2(2,:), 'r')
    plot(LS.M3(2,:), 'k--')
subplot(2, 2, 3)
    hold on
    title('8 steps ahead')
    plot(LS.M1(3,:), 'b')
    plot(LS.M2(3,:), 'r')
    plot(LS.M3(3,:), 'k--')
hold off
legend(strModels)
