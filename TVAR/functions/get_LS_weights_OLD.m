function vecWeights = get_LS_weights_OLD(vecData, strSampleType, strWeightType)
% get_LS_weights
%
% Calculates weights to be used to re-weight a (set of) predictive log
% scores. The weight for time t is only a function of the observation Y(t)
% and the sample on which the model has been estimated [Y(t0), ...,
% Y(t-1)], so the same weights are used for all models.
% 
% Reference: Amisano and Giacomini, Comparing density forecasts via
% weighted likelihood ratio tests, JBES 2007 (AG).
%
% INPUTS:   - vecData: vector of observations
%
%           - strSampleType: string specifying whether the model(s) have
%               been estimated on a recursive or a rolling sample.
%               Acceptable values: 'recursive', 'rolling'.
%           
%           - strWeightType: type of weights to be calculated, i.e. region 
%               of the density on which the evaluation should focus (see AG
%               for details). Acceptable values: 'centre', 'tails',
%               'ltail', 'rtail'.
%
% OUTPUTS:  - vecWeights: vector of weights (for the log-scores)
%
% NOTES As of 19/10, the function works assuming that:
% i) maxfile stays fixed at 56, and
% ii) there is a data0.mat file with the initial (1, ..., T00) sample
% 
% This version: P Alessandri, Oct 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global minfile maxfile T00

% Total number of observations:
TT = length(vecData);

% Output 
vecWeights = NaN(1, TT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get standardised data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vecDataStar = NaN(TT, 1);

for tt = T00+minfile+1 : T00+maxfile

    if strcmpi(strSampleType, 'recursive')    
    vecDataStar(tt) = ( vecData(tt) - mean(vecData(1:tt-1)) ) / std(vecData(1:tt-1));
    
    elseif strcmpi(strSampleType, 'rolling')
    vecDataStar(tt) = ( vecData(tt) - mean(vecData(tt-T00:tt-1)) ) / std(vecData(tt-T00:tt-1));
    
    else
        error('Unknown kind of sample: check inputs')
    end
end

% This is Yst in Amisano-Giacomini equation 1, p.179.
% Note that we standardise Y(t) using moments calculated on data up to t-1. 

% % % Get rid of NaNs:
% % vecDataStar = vecDataStar(T00+1:end);

%% Calculate weights using standard normal distribution:
switch strWeightType
    
    case 'centre'
        vecWeights = pdf('Normal', vecDataStar, 0, 1);
        
    case 'tails'
        vecWeights = 1 - pdf('Normal', vecDataStar, 0, 1) / pdf('Normal', 0, 0, 1) ;
    
    case 'ltail'
        vecWeights = 1 - cdf('Normal', vecDataStar, 0, 1);
        
    case 'rtail'
        vecWeights = cdf('Normal', vecDataStar, 0, 1);
        
    otherwise
        error('Unknown weighting function: check inputs')
end

% Trim vector to remove NaNs:
vecWeights = vecWeights(T00+minfile+1:end);

end
