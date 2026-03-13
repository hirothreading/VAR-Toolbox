%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE over time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'RMSE, horizon 1, all')
subplot(1, 2, 1)
plot(VAR.rmse(:,:,1)), title('VAR'), xlabel('Time') % , ylim([-30 5])
subplot(1, 2, 2)
plot(TAR.rmse(:,:,1)), title('TAR'), xlabel('Time') % , ylim([-30 5])
legend(vnamesvar)

figure('Name', 'RMSE, horizon max, all')
subplot(1, 2, 1)
plot(VAR.rmse(:,:,end)), title('VAR'), xlabel('Time') % , ylim([-20 5])
subplot(1, 2, 2)
plot(TAR.rmse(:,:,end)), title('TAR'), xlabel('Time') % , ylim([-20 5])
legend(vnamesvar)

% 
% figure('Name', 'LS, horizon 1')
% for ii = 1:4
%     subplot(2, 2, ii)
%     hold on
%     plot(VAR.ls(:,ii,1), 'b'), axis square, title(vnamesvar(ii))
%     plot(TAR.ls(:,ii,1), 'r'), axis square, title(vnamestar(ii))
%     hold off
%     xlabel('Time')
% end
% legend('VAR', 'TAR')
% 
% if length(horizons)>=2
% figure('Name', 'LS, horizon 2')
% for ii = 1:4
%     subplot(2, 2, ii)
%     hold on
%     plot(VAR.ls(:,ii,2), 'b'), axis square, title(vnamesvar(ii))
%     plot(TAR.ls(:,ii,2), 'r'), axis square, title(vnamestar(ii))
%     hold off
%     xlabel('Time')
% end
% legend('VAR', 'TAR')
% end
% 
% if length(horizons)>=3
% figure('Name', 'LS, horizon 3')
% for ii = 1:4
%     subplot(2, 2, ii)
%     hold on
%     plot(VAR.ls(:,ii,3), 'b'), axis square, title(vnamesvar(ii))
%     plot(TAR.ls(:,ii,3), 'r'), axis square, title(vnamestar(ii))
%     hold off
%     xlabel('Time')
% end
% legend('VAR', 'TAR')
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% clean up
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% delete temp*
