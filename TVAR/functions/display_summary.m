% Display_summary:
% Script to be run after evaluate_forecasts, to show summary stats on
% screen. P.A. Jan 2013.

% Create strings to label all model comparisons as 'ModelX-ModelY':
strComp_2_1 = strcat(models{2},'-',models{1});
strComp_3_1 = strcat(models{3},'-',models{1});
strComp_3_2 = strcat(models{3},'-',models{2});


disp('------------------------------------------------------------------')
disp('ABSOLUTE VALUES')
disp('------------------------------------------------------------------')

disp('                            RMSE')
printmat(VAR.M.rmse', models{1}, cell2mat(vnamesvar), num2str(horizons));
printmat(TAR.M.rmse', models{2}, cell2mat(vnamestar), num2str(horizons));
printmat(TVTP.M.rmse',models{3}, cell2mat(vnamestar), num2str(horizons));

disp('                            LOG SCORES')
printmat(VAR.M.ls', models{1}, cell2mat(vnamesvar), num2str(horizons));
printmat(TAR.M.ls', models{2}, cell2mat(vnamestar), num2str(horizons));
printmat(TVTP.M.ls',models{3}, cell2mat(vnamestar), num2str(horizons));

% disp('------------------------------------------------------------------')
% disp('Difference in scores')
% disp('------------------------------------------------------------------')
% % NB these are not the test statistics. Not useful - skip. 
% 
% disp('                            RMSE stats')
% printmat(TAR.M.rmse' - VAR.M.rmse', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
% printmat(TVTP.M.rmse'- VAR.M.rmse', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
% printmat(TVTP.M.rmse'- TAR.M.rmse', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));
% 
% disp('                            LS stats')
% printmat(TAR.M.ls' - VAR.M.ls', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
% printmat(TVTP.M.ls'- VAR.M.ls', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
% printmat(TVTP.M.ls'- TAR.M.ls', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));

disp('------------------------------------------------------------------')
disp('UNCONDITIONAL TESTS')
disp('------------------------------------------------------------------')

disp('                            RMSE pvalues')
printmat(GWpv.UC.TAR.rmse', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.UC.TVTP.rmse', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.UC.TVTP2.rmse', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));

disp('                            LS pvalues')
printmat(GWpv.UC.TAR.ls', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.UC.TVTP.ls', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.UC.TVTP2.ls', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));

disp('------------------------------------------------------------------')
disp('CONDITIONAL TESTS')
disp('------------------------------------------------------------------')

disp('                            RMSE pvalues')
printmat(GWpv.C.TAR.rmse', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.C.TVTP.rmse', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.C.TVTP2.rmse', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));

disp('                            LS pvalues')
printmat(GWpv.C.TAR.ls', strComp_2_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.C.TVTP.ls', strComp_3_1, cell2mat(vnamesvar), num2str(horizons));
printmat(GWpv.C.TVTP2.ls', strComp_3_2, cell2mat(vnamesvar), num2str(horizons));


disp('------------------------------------------------------------------')
disp('WEIGHTED LOG SCORES (horizon=1)')
disp('------------------------------------------------------------------')
tempwls = [VAR.M.wlsL ; VAR.M.wlsLR ];
printmat(tempwls, models{1}, 'LeftT BothTs', cell2mat(vnamesvar));

tempwls = [TAR.M.wlsL ; TAR.M.wlsLR ];
printmat(tempwls, models{2}, 'LeftT BothTs', cell2mat(vnamestar));

tempwls = [TVTP.M.wlsL ; TVTP.M.wlsLR ];
printmat(tempwls, models{3}, 'LeftT BothTs', cell2mat(vnamestar));

% tempPvsUncond = [ GWpv.UC.wlsL ; GWpv.UC.wlsLR];
% tempPvsCond   = [ GWpv.C.wlsL ;  GWpv.C.wlsLR];
% printmat(tempPvsUncond, 'p-values (unconditional)', 'LeftT BothTs', cell2mat(vnamesvar));
% printmat(tempPvsCond, 'p-values (conditional)', 'LeftT BothTs', cell2mat(vnamesvar));


clear temp*