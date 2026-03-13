function vecPITs = get_PITs(matPaths, vecData)

global maxfile T00

% Initialise vector of log-scores:
vecPITs = NaN(1, maxfile);

% Loop over time to get log-score for each predictive density
for ii = 1:maxfile-1

    vecPITs(ii) = ksdensity(matPaths(:,ii), vecData(T00 + ii), 'function', 'cdf') ;

end
% Note: the last intHorizon forecasts cannot be scored (the corresponding
% ob