%% Housekeeping
close all
clear
clc
%% Load Data
[DataMat1,XlsInfo]=xlsread('../xls_files/Dataset.xlsx','Dataset_levels');
BeginEstim=find_variable_indices({'31/03/1975'},XlsInfo(2:end,1));
EndEstim=find_variable_indices({'31/12/2008'},XlsInfo(2:end,1));
DataMat2=DataMat1(BeginEstim:EndEstim,5:end);
DataMat2=DataMat2(:,[1 2 4 3 5 6 7 8 9]);
DataMat2(:,1)=100*log(DataMat2(:,1)./100);
DataMat2(:,4)=100*log(DataMat2(:,4));
DataMat3=DataMat2(1:end,1:9);
% DataMat3(:,1)=DataMat2(2:end,1)-DataMat2(1:end-1,1);
% DataMat3(:,3)=DataMat2(2:end,3)-DataMat2(1:end-1,3);
% DataMat3(:,4)=DataMat2(2:end,4)-DataMat2(1:end-1,4);
% DataMat3(:,5)=DataMat2(2:end,5)-DataMat2(1:end-1,5);
% DataMat3(:,6)=DataMat2(2:end,6)-DataMat2(1:end-1,6);

dy=cols(DataMat3);
T=rows(DataMat3);
%% SoE VAR Estimation
PriorBVARMinnesota.burn=1000;
PriorBVARMinnesota.ndraw=2000;
PriorBVARMinnesota.kskip=1;
PriorBVARMinnesota.lamda=2; %Looser (?) prior
PriorBVARMinnesota.muC=0.001;
% PriorBVARMinnesota.tau=10*PriorBVARMinnesota.lamda;
PriorBVARMinnesota.const=1;
PriorBVARMinnesota.nolags=4;
PriorBVARMinnesota.lintrend=1;
PriorBVARMinnesota.BlockIndex=[1 5;4 9];
warning off all
PostBVARMinnesota = bvar_minnesota_estim(PriorBVARMinnesota,DataMat3);
warning on all
%% VAR Identification
SignRes=zeros(4,4);
%| Flow Supply Shock
SignRes(1,1)=-1;
SignRes(2,1)=-1;
SignRes(3,1)=1;


%| Flow Demand Shock
SignRes(1,2)=1;
SignRes(2,2)=1;
SignRes(3,2)=1;

%| Speculative Demand Shock
SignRes(1,3)=1;
SignRes(2,3)=-1;
SignRes(3,3)=1;
SignRes(4,3)=1;
S=size(PostBVARMinnesota.AutoParmChain,3);


SuppIRFs=zeros(dy,12,S);
DemIRFs=zeros(dy,12,S);
SpecIRFs=zeros(dy,12,S);
TermA=0;
for j = 1 : S
    j
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);
    [Rmat,FindIndex]=kulian_jae_identification(BetaMat,SigmaMat,SignRes);
    if FindIndex
        TermA=TermA+1;
        SuppIRFs(:,:,TermA)=varirf(BetaMat,Rmat(:,1),12);
        DemIRFs(:,:,TermA)=varirf(BetaMat,Rmat(:,2),12);
        SpecIRFs(:,:,TermA)=varirf(BetaMat,Rmat(:,3),12);
    end
end

DistSuppIRFs=quantile(SuppIRFs(:,:,1:TermA),[0.5 0.16 0.84],3);
DistDemIRFs=quantile(DemIRFs(:,:,1:TermA),[0.5 0.16 0.84],3);
DistSpecIRFs=quantile(SpecIRFs(:,:,1:TermA),[0.5 0.16 0.84],3);

%% Plots
VarNames={'Oil Production';'Activity';'Oil Price';'Invetnories';...
    'GDP';'CPI';'Policy Rate';'Exchange Rate';'Long Rate'};
xxx=1:12;
figure('name','supply')
for i = 1 : dy
    subplot(3,3,i)
    plotx6(xxx,squeeze(DistSuppIRFs(i,:,:)))
    axis tight
    title(VarNames{i})
end

figure('name','demand')
for i = 1 : dy
    subplot(3,3,i)
    plotx6(xxx,squeeze(DistDemIRFs(i,:,:)))
    axis tight
    title(VarNames{i})
end

figure('name','speculation')
for i = 1 : dy
    subplot(3,3,i)
    plotx6(xxx,squeeze(DistSpecIRFs(i,:,:)))
    axis tight
    title(VarNames{i})
end
