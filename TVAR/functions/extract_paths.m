function matData = extract_paths(strFcFolder, strModel, intVariable, intHorizon)
% -------------------------------------------------------------------------
% extract_draws:
%
% the function extracts draws for a specific variable and horizon from a 
% specific model. 
% 
% INPUTS:   - strFcFolder: string identifying subfolder where the forecasts
%               are stored. That will be {main code folder}\strFcFolder
%           - strModel: string identifying the model
%           - intVariable: integer number identifying the variable
%           - intHorizon: integer number defining the horizon of the
%           forecast
%
% OUTPUTS:  - matData: (REPS)*(No. Samples) matrix of simulated values for 
%               variable number (intVariable) in model (strModel). The
%               values are (intHorizon)-steps ahead forecasts.
%
% P Alessandri, Sept 2012.
% 19/10: changed 'file' loop to take account of minfile=0
% -------------------------------------------------------------------------

global strFolder minfile maxfile REPS BURN HORZ

% Check inputs: -----------------
strAllModels = {'benchmark', 'asset', 'bank'};
if sum(strcmpi(strModel, strAllModels)) == 0
    error('Unknown model: check inputs')
end

if intHorizon > HORZ
    error('Number of steps-ahead exceeds simulation horizon: check inputs')
end
% -------------------------------

% Initialise output:
arrData = NaN(REPS-BURN, HORZ, maxfile-minfile+1);
matData = NaN(REPS-BURN, maxfile-minfile+1);
% NOTE. The dimension is (maxfile-minfile+1) because the labelling 
% convention of the data* and forecast* output files is such that:
%
%   forecast{min}   = E(y/data{min})    = E(y/y0 ... yT00) (assuming min=0)
%   forecast{min+1} = E(y/data{min+1})  = E(y/y0 ... yT00+1)
%	...
%   forecast{Max}   = E(y/data{Max})    = E(y/y0 ... yT00+Max)
%
% So we have Max-Min+1 forecasts. Eg we get a single forecast if Max=Min 
% (in which case only one subsample is used).

% Based on chosen model, define folder where to load paths from:
strOutputFolder = strcat(strFolder, filesep, strFcFolder, filesep, strModel, filesep);

for file=minfile:maxfile
    
    fname=strcat(strOutputFolder,'forecast',num2str(file),'.mat');
    load(fname);
    
    % Extract S*T*1 array of daws for selected variable:
    matData = squeeze(fsave(:, :, intVariable));
    
    % Store data:
    arrData(:, :, file+1) = matData;
    % The +1 shift means that {(:,:,1), (:,:,2), ...} are the paths
    % contained in {forecast0, forecast1, ...}, which are generated using 
    % {data0, data1, ...}.
    
end

% Extract matrix with values for the chosen horizon:
matData = squeeze(arrData(:, intHorizon, :));