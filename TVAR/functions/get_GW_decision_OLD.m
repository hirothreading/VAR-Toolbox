function drule = get_GW_decision_OLD(arrDscores)
% Calculates the "decision rule" suggested by Giacomini-White 2006, p.1558.
% The instruments are assumed to be h = [1 dL(t)]

% Get dimensions:
[T N H] = size(arrDscores);

drule = NaN(T, N, H);

% Size of smallest (initial) sample:
T0 = 12;

% Loop over horizons, variables, time:
for hh = 1:H
    for ii = 1:N
        for tt = T0+1 : T
            
            vecDL = arrDscores(1:tt, ii, hh);
            tempX = [ones(tt-1, 1) vecDL(1:tt-1)];
            tempY = vecDL(2:end);

            beta = tempX\tempY;
            tempFitted = beta' * tempX';


            drule(tt, ii, hh) = tempFitted(end);
            
        end     
    end
end


end

%% TOPPATO: regression must be estimated at every t, so that coefficients 
% only incorporate time-t information. Here we are getting a single
% estimate of delta which incorporates full-sample information.