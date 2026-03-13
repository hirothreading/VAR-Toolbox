%% Housekeeping
close all
clear
clc
%% Load Data
[DataMat1,XlsInfo]=xlsread('../xls_files/Dataset3.xlsx','Dataset_levels');
BeginEstim=find_variable_indices({'31/03/1986'},XlsInfo(3:end,1));
EndEstim=find_variable_indices({'31/12/2015'},XlsInfo(3:end,1));
DataMat2=DataMat1(BeginEstim:EndEstim,:);
VarNames=XlsInfo(2,2:end);
VarNames=VarNames(:,[9 4 17 8 10 16 11 12 13 18]);
DataMat3=DataMat2(:,[9 4 17 8 10 16 11 12 13 18]);
DataMat=DataMat3(1:end,:);
DataMat(:,1)=cumsum(DataMat(:,1));
% DataMat(:,2:9)=DataMat3(2:end,2:9)-DataMat3(1:end-1,2:9);

dy=cols(DataMat);
T=rows(DataMat);
BayesianEstim=1;
NoIrfs=40;
NoHorizon=4;
%% VAR
VAR.p=1; % Number of Lags
VAR.irhor=40; % Impulse Response Horizon
Y2=DataMat(:,2:end);
VAR.vars=Y2; %order the first variable in data is where the shock identified
VAR.proxies = DataMat(:,1); %proxy variable

VAR01=VAR;
VAR01=doProxySVAR_single(VAR01);
temp=VAR01.irs;
ProxyIRFs=temp';

VARbs = doProxySVARbootstrap_single(VAR01,1000,90);
%% Plots
xxx=1:NoIrfs;
figure('name','Supply-IRFs')
for i = 1 : dy-1
    subplot(3,4,i)
    plotx1(xxx(1:12),[ProxyIRFs(i,1:12)' VARbs.irsL(1:12,i) VARbs.irsH(1:12,i)])
    axis tight
    title(VarNames{i+1})
end
