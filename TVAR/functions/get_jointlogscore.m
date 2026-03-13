function vecLogScores = get_jointlogscore(matPathsX, matPathsY, vecDataX, vecDataY, vecHorizons)
% -------------------------------------------------------------------------
% get_jointlogscore:
%
% Calculates join log-score for for a set of predicted values and actual
% observations on two variables (X, Y)
% 
% INPUTS:   - matPathsX: (# draws)*(# periods) matrix of simulated values for
%               variable X.
%           - matPathsY: (# draws)*(# periods) matrix of simulated values for
%               variable Y.
%           - vecDataX: (# periods)*1 matrix of observations on X
%           - vecDataY: (# periods)*1 matrix of observations on Y
%           - vecHorizons: vector defining the horizons at which the
%               forecasts are evaluated
%           - strType: 'plain' to score period-specific forecasts, 'cumul' to
%               score cumulative forecasts. Choice depends on the variable:
%               use 'cum' for changes/growth rates.
%
% OUTPUTS:  - vecLogScores: 1*... vector of log-scores. The value
%               is NaN for the densities for which we do not have
%               observations (ie the last ones).
%
% P Alessandri, Jan 2014. Note: this version (unlike the 7 Jan 2014 one)
% does not need to discriminate between "norm" and "cumul" cases.
% -------------------------------------------------------------------------


% Initialise vector of log-scores:
vecLogScores = NaN(1, length(vecHorizons));

% Loop over horizons to get log-score for each predictive density

for i = 1:length(vecHorizons)
    
    ii=vecHorizons(i);
    
    pathsX = matPathsX(:,ii);
    pathsY = matPathsY(:,ii);
    obsX   = vecDataX(ii);
    obsY   = vecDataY(ii);            
    
    % Estimate joint pdf on a grid: 
    [bandwidth, pdf, X, Y] = kde2d([pathsX pathsY], 100);

	% -------------------------------------------------------------------
	% Replace tiny negatives (due to numerical issues) with eps
    if min(min(pdf))<-0.00001
        warning('kde2d returns negative pdf')
    end
    pdf_ind = (pdf<0);
    pdf(pdf_ind) = eps;
	% -------------------------------------------------------------------
    
    % Get the value of the pdf in (obsX, obsY) by interpolation:
    pdf0 = interp2(X, Y, pdf, obsX, obsY);
    
	% -------------------------------------------------------------------
	% Replace NaN due to the observation being outside the domain of the
	% pdf with eps:
    if isnan(pdf0)
        if obsX<min(min(X)) | obsX>max(max(X)) | obsY<min(min(Y)) | obsY>max(max(Y)) 
            pdf0 = eps;
        else
            error('pdf interpolation returns NaN')
        end
    end
	% -------------------------------------------------------------------
    
    
    
    vecLogScores(i) = log(pdf0);
    
    end

end