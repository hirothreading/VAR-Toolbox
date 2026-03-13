function matData = get_data(strDataFolder, strModel)
% -------------------------------------------------------------------------
% get_data:
%
% ... 
% 
% INPUTS:   - 
%
% OUTPUTS:  - 
%
% P Alessandri, Sept 2012
% -------------------------------------------------------------------------

global strFolder minfile maxfile

% Check inputs:
strAllModels = {'benchmark', 'asset', 'bank'};
if sum(strcmpi(strModel, strAllModels)) == 0
    error('Unknown model: check inputs')
end

% Initialise output:
% matData = NaN;

% Based on chosen model, define folder where to load data from:
strTemp = strcat(strDataFolder, strModel);
strFullFolderPath = strcat(strFolder, filesep, strDataFolder, filesep, strTemp, filesep);
% 1st line needed to merge 'data' and 'model' (eg \databenchmark\)

% Define name of 'data' file to be loaded, and load it:
fname = strcat(strFullFolderPath,'data',num2str(maxfile),'.mat');
load(fname);
% We pick the last one, which has all available observations.

% Define output:
matData = dataout;

