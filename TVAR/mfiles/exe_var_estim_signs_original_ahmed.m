%% Housekeeping
%close all
clear
clc
%% Load Data
Country = 'US'; %enter the country of interest; make sure the excel file exists in the xls_files folder

%%Because of different spreadsheet names for the UK and other Countries
if strcmpi(Country, 'UK')%Country == 'UK'
    spreadsheet = 'Dataset_levels';
else
    spreadsheet = 'Transformed';
end

%%loads data for the relevant country over the relevant time period
[DataMat1,XlsInfo]=xlsread(strcat('../xls_files/',Country,'.xlsx'),spreadsheet);
BeginEstim=find_variable_indices({'31/03/1986'},XlsInfo(3:end,1)); %1986
EndEstim=find_variable_indices({'31/12/2015'},XlsInfo(3:end,1));
DataMat2=DataMat1(BeginEstim:EndEstim,:);
VarNames=XlsInfo(2,2:end);

%%selects macro-variables of interest; to be chosen for each country;
%exclude variables with missing data (see each excel file individually)
if strcmpi(Country, 'UK')%Country == 'UK'
    VarNames=VarNames(:,[9 4 19 8   10 16 11 12 13 14 18]);   %17=OECD IP; 19=World IP; 5=REA; 9=Surprise
    DataMat3=DataMat2(:,[9 4 19 8   10 16 11 12 13 14 18]);
else
    VarNames=VarNames(:,[6 1 14 5   8  10 11 12 13]); %6=Surpr.; 2=Activity Index; 7=OECD IP; 14=World IP; 3=EIA Invent.
    DataMat3=DataMat2(:,[6 1 14 5   8  10 11 12 13]);
end

DataMat=DataMat3(1:end,:);

%%uncomment if you include data on energy surprises
DataMat(:,1)=cumsum(DataMat(:,1)); 

%%uncomment if you want to run the model with first-different
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
PriorBVARMinnesota.BlockIndex=[1 5;4 9];
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

SuppIRFs=zeros(dy+1,NoIrfs,S);
SuppFVDs=zeros(dy+1,NoIrfs,S);

DemIRFs=zeros(dy+1,NoIrfs,S);
DemFVDs=zeros(dy+1,NoIrfs,S);

OilDemIRFs=zeros(dy+1,NoIrfs,S);
OilDemFVDs=zeros(dy+1,NoIrfs,S);

Sdmat=zeros(2,dy*PriorBVARMinnesota.nolags);
Sdmat(1,5)=1;   %1stLag vector; output
Sdmat(2,6)=1;   %1stLag vector; cpi
Sdmat(2,33)=-1; %4thLag vector; cpi
    
RmatStore=zeros(dy,dy,S);
normval=zeros(S,1);

yHD=zeros(rows(DataMat)-4,dy+1);
cpiHD=zeros(rows(DataMat)-4,dy+1);

[Y,X]=makelagsmat(DataMat,PriorBVARMinnesota.nolags);
X=[ones(rows(Y),1) [0:rows(Y)-1]' X];


for j = 1 : S
    j
    
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);
    [Rmat,FindIndex,normval(j)]=kulian_jae_identification_original(SigmaMat,signrest(:));    %
    RmatStore(:,:,j)=Rmat;
    
    SuppIRFs(1:dy,:,j)=varirf(BetaMat,Rmat(:,1),NoIrfs);
    DemIRFs(1:dy,:,j)=varirf(BetaMat,Rmat(:,2),NoIrfs);
    OilDemIRFs(1:dy,:,j)=varirf(BetaMat,Rmat(:,3),NoIrfs);
    
    [~,Dmat] = forecast_var_decomposition(BetaMat,Rmat,[1:NoIrfs]');
    
    SuppFVDs(1:dy,:,j)=squeeze(Dmat(1:end-1,1,:))';
    DemFVDs(1:dy,:,j)=squeeze(Dmat(1:end-1,2,:))';
    OilDemFVDs(1:dy,:,j)=squeeze(Dmat(1:end-1,3,:))';
    
    
    bigpi=varcompanion(BetaMat);
    bigw=zeros(dy*PriorBVARMinnesota.nolags,dy);
    bigw(1:dy,:)=Rmat;
    
    obsirf  = dsge_obs_irf(bigpi,bigw,Sdmat,NoIrfs);
    [~,Dmat] = forecast_dsge_decomposition(Sdmat,bigpi,bigw,[1:NoIrfs]');
    
    SuppIRFs(dy+1,:,j)  =obsirf(2,:,1);
    DemIRFs(dy+1,:,j)   =obsirf(2,:,2);
    OilDemIRFs(dy+1,:,j)=obsirf(2,:,3);
    
    SuppFVDs(dy+1,:,j)  =Dmat(1:end-1,1,2)';
    DemFVDs(dy+1,:,j)   =Dmat(1:end-1,2,2)';
    OilDemFVDs(dy+1,:,j)=Dmat(1:end-1,3,2)';
    
    
    errormat=inv(Rmat)*(Y'-PostBVARMinnesota.AutoParmChain(:,:,j)*X');
    
    DecompVAR=var_historical_decomp(PostBVARMinnesota.AutoParmChain(:,:,j),...
        Rmat,DataMat,errormat,1,1);
    
    yHD=yHD+[DecompVAR.yfs(:,:,5) DecompVAR.initeffecty(5,:)'];
    cpiHD=cpiHD+[DecompVAR.yfs(:,:,6) DecompVAR.initeffecty(6,:)'];
    
end

DistSuppIRFs=quantile(SuppIRFs,[0.5 0.16 0.84],3);
DistSuppFVDs=quantile(SuppFVDs,[0.5 0.16 0.84],3);

DistDemIRFs=quantile(DemIRFs,[0.5 0.16 0.84],3);
DistDemFVDs=quantile(DemFVDs,[0.5 0.16 0.84],3);

DistOilDemIRFs=quantile(OilDemIRFs,[0.5 0.16 0.84],3);
DistOilDemFVDs=quantile(OilDemFVDs,[0.5 0.16 0.84],3);
yHD=yHD./S;
cpiHD=cpiHD./S;

%% Kilian's Unique Matrix

% [~,minIndex]=min(normval);
% RmatUnique=RmatStore(:,:,minIndex);
% BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,minIndex);
%
% SuppIRFsUnique=varirf(BetaMat,RmatUnique(:,1),NoIrfs);
% DemIRFsUnique=varirf(BetaMat,RmatUnique(:,2),NoIrfs);
% OilDemIRFsUnique=varirf(BetaMat,RmatUnique(:,3),NoIrfs);
%
% [~,Dmat] = forecast_var_decomposition(BetaMat,RmatUnique,[1:NoIrfs]');
%
% SuppFVDsUnique=squeeze(Dmat(1:end-1,1,:))';
% DemFVDsUnique=squeeze(Dmat(1:end-1,2,:))';
% OilDemFVDsUnique=squeeze(Dmat(1:end-1,3,:))';
%
%
% bigpi=varcompanion(BetaMat);
% bigw=zeros(dy*PriorBVARMinnesota.nolags,dy);
% bigw(1:dy,:)=RmatUnique;
%
% obsirf  = dsge_obs_irf(bigpi,bigw,Sdmat,NoIrfs);
% [~,Dmat] = forecast_dsge_decomposition(Sdmat,bigpi,bigw,[1:NoIrfs]');
%
% SuppIRFsUnique(dy+1,:)  =obsirf(2,:,1);
% DemIRFsUnique(dy+1,:)   =obsirf(2,:,2);
% OilDemIRFsUnique(dy+1,:)=obsirf(2,:,3);
%
% SuppFVDsUnique(dy+1,:)  =Dmat(1:end-1,1,2)';
% DemFVDsUnique(dy+1,:)   =Dmat(1:end-1,2,2)';
% OilDemFVDsUnique(dy+1,:)=Dmat(1:end-1,3,2)';
%
% [Y,X]=makelagsmat(DataMat,PriorBVARMinnesota.nolags);
% X=[ones(rows(Y),1) [0:rows(Y)-1]' X];
%
% errormat=inv(RmatUnique)*(Y'-PostBVARMinnesota.AutoParmChain(:,:,minIndex)*X');
% DecompVAR=var_historical_decomp(PostBVARMinnesota.AutoParmChain(:,:,minIndex),...
%     Rmat,DataMat,errormat,1,1);
%

%% Plots
VarNames=[VarNames {'Inflation'}];
xxx=1:NoIrfs;
figure('name','Supply-IRFs')
for i = 1 : dy+1
    subplot(3,4,i)
    plotx6(xxx(1:20),squeeze(DistSuppIRFs(i,1:20,:)))
%    hold on
%    plot(xxx(1:20),SuppIRFsUnique(i,1:20),'-.r','LineWidth',4)
%    hold off
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
for i = 1 : dy+1
    subplot(3,4,i)
    plotx6(xxx(1:20),squeeze(DistDemIRFs(i,1:20,:)))
%     hold on
%     plot(xxx(1:20),DemIRFsUnique(i,1:20),'-.r','LineWidth',4)
%     hold off
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

figure('name','Oil-Demand-IRFs')
for i = 1 : dy
    subplot(3,4,i)
    plotx6(xxx(1:20),squeeze(DistOilDemIRFs(i,1:20,:)))
%     hold on
%     plot(xxx(1:20),OilDemIRFsUnique(i,1:20),'-.r','LineWidth',4)
%     hold off
    axis tight
    title(VarNames{i})
end
%
% figure('name','Oil-Demand-FVDs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx,squeeze(DistOilDemFVDs(i,:,:)))
%     axis tight
%     title(VarNames{i})
% end

%% xls
% xlswrite('../xls_files/HDs2_US.xlsx',[yHD DataMat(5:end,5)],'yHD','B58:L173');
% 
% xlswrite('../xls_files/HDs2_US.xlsx',[cpiHD DataMat(5:end,6)],'cpiHD','B58:L173');