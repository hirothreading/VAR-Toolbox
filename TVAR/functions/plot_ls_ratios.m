%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIFFERENCES between  log-scores over time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

???     ???     ??? 
??? ELIMINATE   ??? 
???     ??? 	??? 

%% Check if LSs are all negative
tempCheckY  = sum(sum(sum(VAR.ls(:,1,:)>0))) + sum(sum(sum(TAR.ls(:,1,:)>0))) + sum(sum(sum(TVTP.ls(:,1,:)>0))) ;
tempCheckPi = sum(sum(sum(VAR.ls(:,2,:)>0))) + sum(sum(sum(TAR.ls(:,2,:)>0))) + sum(sum(sum(TVTP.ls(:,2,:)>0))) ;
tempCheckR  = sum(sum(sum(VAR.ls(:,3,:)>0))) + sum(sum(sum(TAR.ls(:,3,:)>0))) + sum(sum(sum(TVTP.ls(:,3,:)>0))) ;
tempCheckF  = sum(sum(sum(VAR.ls(:,4,:)>0))) + sum(sum(sum(TAR.ls(:,4,:)>0))) + sum(sum(sum(TVTP.ls(:,4,:)>0))) ;
tempCheckYPi= sum(sum(sum(VAR.jls(:,:)>0)))  + sum(sum(sum(TAR.jls(:,:)>0))) + sum(sum(sum(TVTP.jls(:,:)>0))) ;

disp(sprintf('Check for log-scores>0:   Y %s    Pi %s   R %s    F %s    Y&Pi %s', ...
            num2str(tempCheckY), num2str(tempCheckPi), num2str(tempCheckR), num2str(tempCheckF), num2str(tempCheckYPi )))
% ...

%% Other preliminaries 

% Tickmarks for dates:
datesind = 2:48:342;
% Define data and un-do IP and CPI transformation (see arrangedata):
tempData1   = datavar;
tempData1   = tempData1 ./ repmat([12 1 12 1], T,1);
% Get cumulative data
tempData2 = tempData1;
for tt = T0+hh:size(tempData1, 1) 
    tempData2(tt, [1 3])  = sum(tempData1(tt-hh+1:tt, [1 3]), 1);
end


%% All variables individually, all horizons, VAR vs TAR

figure('Name', 'LS diffs, all vars, TAR vs VAR')

tempLineW = 1.5;

for ii=1:4
    
    subplot(2,2,ii)
    tempX = squeeze(TAR.ls(:,ii,:) - VAR.ls(:,ii,:));
    hold on
%     plot(datesn2, tempX(:,1), 'g', 'LineWidth', tempLineW);
    plot(datesn2, tempX(:,2), 'b', 'LineWidth', tempLineW);
%     plot(datesn2, tempX(:,3), 'c', 'LineWidth', tempLineW);
    plot(datesn2, tempX(:,4), 'r', 'LineWidth', tempLineW);
    hline(0, 'k')
    hold off
    set(gca, 'XTick', datesn2(datesind));
    set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
    axis tight
    set(gca, 'ylim', [median(tempX(:,1))-2*std(tempX(:,1)) median(tempX(:,1))+2*std(tempX(:,1))])
    title(vnamesvar(ii))

end

% legend('1M', '3M', '6M', '12M', 'Location', 'Best'), legend BOXOFF
legend('3M', '12M', 'Location', 'Best'), legend BOXOFF

%% Joint (Y, Pi), all horizons, VAR vs TAR
figure('Name', 'LS diffs, (y,pi), TAR vs VAR')

tempX = squeeze(TAR.jls(:,:) - VAR.jls(:,:));

hold on
plot(datesn2, tempX(:,1), 'g', 'LineWidth', tempLineW);
plot(datesn2, tempX(:,2), 'b', 'LineWidth', tempLineW);
plot(datesn2, tempX(:,3), 'c', 'LineWidth', tempLineW);
plot(datesn2, tempX(:,4), 'r', 'LineWidth', tempLineW);
hline(0, 'k')
hold off
set(gca, 'XTick', datesn2(datesind));
set(gca, 'XTickLabel', datestr(datesn2(datesind), 'mm/yy'));
axis 'tight'
legend('1M', '3M', '6M', '12M', 'Location', 'NorthWest'), legend BOXOFF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delete temp*
