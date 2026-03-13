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
% DataMat(:,2:8)=DataMat3(2:end,2:8)-DataMat3(1:end-1,2:8);

%dy=cols(DataMat);
%T=rows(DataMat);
dy=size(DataMat,2);
T=size(DataMat,1);
BayesianEstim=1;
NoIrfs=40;
NoHorizon=4;
%% SoE VAR Estimation
PriorBVARMinnesota.burn=1000;
PriorBVARMinnesota.ndraw=2000;
PriorBVARMinnesota.kskip=1;
PriorBVARMinnesota.lamda=2; %Looser (?) prior
PriorBVARMinnesota.muC=0.0001;
% PriorBVARMinnesota.tau=10*PriorBVARMinnesota.lamda;
PriorBVARMinnesota.const=1;
PriorBVARMinnesota.nolags=4;
PriorBVARMinnesota.lintrend=1;
PriorBVARMinnesota.BlockIndex=[1 5;4 10];
warning off all
PostBVARMinnesota = bvar_minnesota_estim(PriorBVARMinnesota,DataMat);
warning on all
S=size(PostBVARMinnesota.AutoParmChain,3);

%% Shock Identification
signrest=zeros(4,4);
signrest(1,1)=1;
signrest(2,1)=-1;
signrest(3,1)=-1;
signrest(4,1)=1;

signrest(1,2)=1;
signrest(2,2)=1;
signrest(3,2)=1;
signrest(4,2)=1;

signrest(1,3)=1;
signrest(2,3)=1;
signrest(3,3)=-1;
signrest(4,3)=1;

SuppIRFs=zeros(dy,NoIrfs,S);
SuppFVDs=zeros(dy,NoIrfs,S);

DemIRFs=zeros(dy,NoIrfs,S);
DemFVDs=zeros(dy,NoIrfs,S);

OilDemIRFs=zeros(dy,NoIrfs,S);
OilDemFVDs=zeros(dy,NoIrfs,S);

for j = 1 : S
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);
    [Rmat,FindIndex]=kulian_jae_identification_original(SigmaMat,signrest(:));    %
    
    SuppIRFs(:,:,j)=varirf(BetaMat,Rmat(:,1),NoIrfs);
    DemIRFs(:,:,j)=varirf(BetaMat,Rmat(:,2),NoIrfs);
    OilDemIRFs(:,:,j)=varirf(BetaMat,Rmat(:,3),NoIrfs);
    
    SuppFVDs(:,:,j)=forecast_var_decomp_noP0mat(SuppIRFs(:,:,j),BetaMat,SigmaMat);
    DemFVDs(:,:,j)=forecast_var_decomp_noP0mat(DemIRFs(:,:,j),BetaMat,SigmaMat);
    OilDemFVDs(:,:,j)=forecast_var_decomp_noP0mat(OilDemIRFs(:,:,j),BetaMat,SigmaMat);
    j
end

DistSuppIRFs=quantile(SuppIRFs,[0.5 0.16 0.84],3);
DistSuppFVDs=quantile(SuppFVDs,[0.5 0.16 0.84],3);

DistDemIRFs=quantile(DemIRFs,[0.5 0.16 0.84],3);
DistDemFVDs=quantile(DemFVDs,[0.5 0.16 0.84],3);

DistOilDemIRFs=quantile(OilDemIRFs,[0.5 0.16 0.84],3);
DistOilDemFVDs=quantile(OilDemFVDs,[0.5 0.16 0.84],3);

%% Plots
xxx=1:NoIrfs;
figure('name','Supply-IRFs')
for i = 1 : dy
    subplot(3,4,i)
    plotx6(xxx(1:20),squeeze(DistSuppIRFs(i,1:20,:)))
    axis tight
    title(VarNames{i})
end

% figure('name','Supply-FVDs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx,squeeze(DistSuppFVDs(i,:,:)))
%     axis tight
%     title(VarNames{i})
% end

figure('name','Demand-IRFs')
for i = 1 : dy
    subplot(3,4,i)
    plotx6(xxx(1:20),squeeze(DistDemIRFs(i,1:20,:)))
    axis tight
    title(VarNames{i})
end

% figure('name','Demand-FVDs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx,squeeze(DistDemFVDs(i,:,:)))
%     axis tight
%     title(VarNames{i})
% end
% 
% figure('name','Oil-Demand-IRFs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx(1:20),squeeze(DistOilDemIRFs(i,1:20,:)))
%     axis tight
%     title(VarNames{i})
% end
% 
% figure('name','Oil-Demand-FVDs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx,squeeze(DistOilDemFVDs(i,:,:)))
%     axis tight
%     title(VarNames{i})
% end