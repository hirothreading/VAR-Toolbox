%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% G-W decision rules 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datesind = 2:48:342;

%% Joint LS decision rule for (Y, PI) at all horizons:

figure('Name', 'Decision rule, JLS')

tempLS = squeeze(GWdr.jls0(:,:)); 

hold on
plot(datesn2, tempLS(:,2), 'b--', 'Linewidth', 1.5)
plot(datesn2, tempLS(:,4), 'r', 'Linewidth', 1.5)
hline(0, 'k')
hold off

axis tight
set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
legend('3M', '12M', 'Location', 'NorthWest'), legend BOXOFF

%% housekeeping
delete temp*

