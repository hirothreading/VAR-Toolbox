% Display_summary:
% Script to be run after evaluate_forecasts, to show summary stats on
% screen. P.A. Jan 2013.

disp('------------------------------------------------------------------')
disp('RMSE')
disp('------------------------------------------------------------------')
printmat(VAR.M.rmse, 'VAR', num2str(horizons), cell2mat(vnamesvar));
printmat(TAR.M.rmse, 'TAR', num2str(horizons), cell2mat(vnamestar));
printmat(basic.M.rmse, 'basic', num2str(horizons), cell2mat(vnamesbasic));


disp('------------------------------------------------------------------')
disp('LOG SCORES')
disp('------------------------------------------------------------------')

printmat(VAR.M.ls, 'VAR', num2str(horizons), cell2mat(vnamesvar));
printmat(TAR.M.ls, 'TAR', num2str(horizons), cell2mat(vnamestar));
printmat(basic.M.ls, 'basic', num2str(horizons), cell2mat(vnamesbasic));


disp('------------------------------------------------------------------')
disp('WEIGHTED LOG SCORES (horizon=1)')
disp('------------------------------------------------------------------')
tempwlsL = [VAR.M.wlsL ; TAR.M.wlsL;[basic.M.wlsL nan(1,1)]];
printmat(tempwlsL, 'Left Tail', 'VAR TAR BASIC', cell2mat(vnamesvar));

tempwlsLR = [VAR.M.wlsLR ; TAR.M.wlsLR; [basic.M.wlsLR nan(1,1)]];
printmat(tempwlsLR, 'Both Tails', 'VAR TAR BASIC', cell2mat(vnamesvar));

clear temp*