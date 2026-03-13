function vecLogScores = get_logscore(matPaths, vecData, intHorizon)
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
% P Alessandri, Oct 2012
%   - 24 October: changed timing convention for the output vector
% -------------------------------------------------------------------------

global minfile maxfile T00

% Initialise vector of log-scores:
vecLogScores = NaN(1, maxfile-minfile);
% NOTE: see extract_paths. Here we use (max-min) instead of (max-min+1)
% because the last set of forecasts E(y/data{Max})= E(y/y0 ... yT00+Max) 
% cannot be evaluated.

% Loop over time to get log-score for each predictive density
for ii = 1 : (maxfile-minfile)-intHorizon+1

    vecLogScores(intHorizon+ii-1) = log( ksdensity(matPaths(:,ii), vecData(T00+intHorizon+ii-1)) );

end
% NOTE:  say horizon=1 and paths(:,t)=[Y(t)/I(t-1)] with t=T00, T00+1, ... Then:
%
%   ii=1        -->  paths(:,1) = Y(T00+1)/I(T00)    is compared to y(T00+1) and stored in vecLS(1)
%   ii=2        -->  paths(:,2) = Y(T00+2)/I(T00+1)  is compared to y(T00+2) and stored in vecLS(2)
% 	...
%   ii=Max-Min-H+1
%               -->  paths(:,d) = Y(T00+d)/I(T00+d-1) is compared to y(T00+d) and stored in vecLS(Max-Min)
% 
% where I(t) is the information set on which the forecast is based and 
% d = Max-Min. It follows that the right matching/comparison is :
%
%       {vecLogScore(t), matPaths(:,t), vecData(T00+t)}
%
% where (t) is the period the forecast refers to (as in matPaths)

end