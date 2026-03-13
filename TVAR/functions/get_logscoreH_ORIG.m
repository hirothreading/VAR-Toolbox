function vecLogScores = get_logscoreH(matPaths, vecData, intHorizon)
% -------------------------------------------------------------------------
% get_logscore:
%
% Calculates log-scores for a set of predicted values and actual
% observations on a specific variable
% 
% INPUTS:   - matPaths: (REPS)*(No. Samples) matrix of simulated values for
%               the variable.
%           - vecData: T*1 vector of actual observations (T>No.Samples)
%           - intHorizon: integer number defining the horizon of the
%               forecast
%
% OUTPUTS:  - vecLogScores: 1*(No.Samples) vector of log-scores. The value
%               is NaN for the densities for which we do not have
%               observations (ie the last ones).
%
% P Alessandri, Sept 2012
% -------------------------------------------------------------------------


% Initialise vector of log-scores:
vecLogScores = NaN(1, length(intHorizon));

% Loop over time to get log-score for each predictive density
for i = 1:length(intHorizon)
    ii=intHorizon(i);
    vecLogScores(i) = log( ksdensity(matPaths(:,ii), vecData(ii)) );

end
% Note: the last intHorizon forecasts cannot be scored (the corresponding
% observations are not available).

end