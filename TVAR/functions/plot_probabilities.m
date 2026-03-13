%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Probabilties of output declines against data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reminder: dimensions are VAR.probs(time, horizon, size of output fall)

% Pick horizon (in terms of months):
hh = 12;

% Choose a size (ordered as in 'probabilities'):
tempsz = 4;

% Tickmarks for dates:
% datesind = 1:12:nobsGR;
datesind = 2:48:342;

% NOTHING TO BE CHANGED FROM HERE ON.

% Get position of the horizon in the results data
hhp = find(horizons==hh);
if isempty(hhp)
    error('There is no such horizon')
end


% Define data and un-do IP and CPI transformation (see arrangedata):
% ------------------------------------------------------------------
tempData1 = datavar;
% Cumulate:
tempData2 = tempData1;
for tt = T0+hh:size(tempData1, 1) 
    tempData2(tt, [1 3])  = sum(tempData1(tt-hh+1:tt, [1 3]), 1);
end
% Divide by horizon to get % per year:
tempData2(:, [1 3]) = tempData2(:, [1 3]) / hh;
% ------------------------------------------------------------------

figure('Name', 'Probs & data ')

% pick percentile:
tempM1 = squeeze(VAR.probs(:, hhp, tempsz ));    
tempM2 = squeeze(TAR.probs(:, hhp, tempsz ));
tempM3 = squeeze(TVTP.probs(:, hhp, tempsz ));

hold on
% area(datesn2, tempData2(T0+hh:T0+T2+hh-1, 1), 'FaceColor',[.7 .8 .9], 'EdgeColor', [1 1 1])

plot(datesn2, tempM3, 'b--', 'LineWidth', 2), axis tight, 
plot(datesn2, tempM1, 'b', 'LineWidth', 2), axis tight,
plot(datesn2, tempM2, 'r', 'LineWidth', 2), axis tight, 

hline(0, 'k')
hold off

xlim ([min(datesn2(datesind)) max(datesn2(datesind))])
set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));

ylim auto
    
% Two percentiles:
legend('VAR^§', 'VAR', 'TAR', 'Location', 'Best'), legend BOXON


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear temp*
