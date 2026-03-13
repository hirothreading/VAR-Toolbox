%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-W decision rules 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datesind = 2:48:342;

hh = 4;

%% LS and RMSE decision rules for a given variable and horizon:
figure('Name', 'Decision rule')

ii = 1; 
tempLS   = squeeze(GWdr.ls(:,ii,hh)); 
tempRMSE = squeeze(GWdr.rmse(:,ii,hh));

hold on
plot(datesn2, tempRMSE, 'b', 'Linewidth', 1.5)
plot(datesn2, tempLS, 'r', 'Linewidth', 1.5)
hline(0, 'k')
hold off

axis tight

% ylim([-1 5])
legend('RMSE', 'LS', 'Location', 'NorthWest'), legend BOXOFF
% title(vnamemsvar(ii))

set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));

%% LS and RMSE decision rules for (IP, CPI):
figure('Name', 'Decision rule (IP & CPI)')

ii = 1; 
tempLS   = squeeze(GWdr.jls(:,hh)); 

hold on
plot(datesn2, tempLS, 'r', 'Linewidth', 1.5)
hline(0, 'k')
hold off

axis tight

% ylim([-1 5])
legend('Joint LS', 'Location', 'NorthWest'), legend BOXOFF
% title(vnamemsvar(ii))

set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));


%% housekeeping
delete temp*

