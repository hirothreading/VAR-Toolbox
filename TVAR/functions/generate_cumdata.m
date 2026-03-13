function matCumData = generate_cumdata(matData, intHorizon) 
% -------------------------------------------------------------------------
% generate_cumdata:
%
% Adds up the rows in the input matrix matData. Supposed to generate
% cumulative growth from an input of period-by-period growth 
% 
% P Alessandri, Oct 2012
% -------------------------------------------------------------------------

% Initialise:
[intS intT] = size(matData);
matCumData = NaN(intS, intT);

% Loop over periods:
for tt = intHorizon : intT
    
    matCumData(:, tt) = sum(matData(:, tt-intHorizon+1:tt), 2);
    
end

% Row (t) contains the cumulative value y(t-h+1)+...+y(t). Eg with horizon
% = 4, cumdata(t) = data(t-3)+data(t-2)+data(t-1)+data(t). 

end
