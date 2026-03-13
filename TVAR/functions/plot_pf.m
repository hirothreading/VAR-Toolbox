%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Point forecasts against data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name', 'PF & data, horizon 1 (all)')

for ii = 1:N
   
    subplot(2, 2, ii)
    hold on
    plot(VAR.pf(:,ii,1), 'b-', 'LineWidth', 2), axis tight,
    plot(TAR.pf(:,ii,1), 'r-', 'LineWidth', 2), axis tight, 
    plot(dataout(T0+1:T0+T2, ii), 'k', 'LineWidth', 1)
    hold off
    title(vnamestar(ii))
    xlabel('Time')
       
end
legend('VAR', 'TAR', 'data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear temp*
