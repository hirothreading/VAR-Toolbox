%% Housekeeping
%close all
clear
clc

%% Load Data
Country = 'Canada'; %enter the country of interest; make sure the excel file exists in the xls_files folder

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
    VarNames=VarNames(:,[9 4 5 8   10  16 11 12 13 14 18]); %1=US_GDP; 5=Activity Index; 17=OECD_IP; 16=UK_IP; 6=EIA Invent.; 9=Surprise
    DataMat3=DataMat2(:,[9 4 5 8   10  16 11 12 13 14 18]);
else
    VarNames=VarNames(:,[6 1 14 5   8  10 11 12 13]); %6=Surpr.; 2=Activity Index; 7=OECD IP; 3=EIA Invent.; 14=World IP
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

SuppIRFs=zeros(dy,NoIrfs,S);
SuppFVDs=zeros(dy,NoIrfs,S);

DemIRFs=zeros(dy,NoIrfs,S);
DemFVDs=zeros(dy,NoIrfs,S);

OilDemIRFs=zeros(dy,NoIrfs,S);
OilDemFVDs=zeros(dy,NoIrfs,S);

for j = 1 : S
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);
    [Rmat,FindIndex]=kulian_jae_identification(SigmaMat,signrest(:));    %
    
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
    subplot(3,3,i)
    plotx6(xxx(1:20),squeeze(DistSuppIRFs(i,1:20,:)))
    axis tight
    title(VarNames{i})
end

% figure('name','Supply-FVDs')
% for i = 1 : dy
%     subplot(3,3,i)
%     plotx6(xxx(1:20),squeeze(DistSuppFVDs(i,1:20,:)))
%     axis tight
%     title(VarNames{i})
% end
% 
figure('name','Demand-IRFs')
for i = 1 : dy
    subplot(3,3,i)
    plotx6(xxx(1:20),squeeze(DistDemIRFs(i,1:20,:)))
    axis tight
    title(VarNames{i})
end

% figure('name','Demand-FVDs')
% for i = 1 : dy
%     subplot(3,3,i)
%     plotx6(xxx(1:20),squeeze(DistDemFVDs(i,1:20,:)))
%     axis tight
%     title(VarNames{i})
% end
% 
% figure('name','Oil-Demand-IRFs')
% for i = 1 : dy
%     subplot(3,4,i)
%     plotx6(xxx(1:40),squeeze(DistOilDemIRFs(i,1:40,:)))
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


%% Concatenating the Country-specific IRFs - to be manually pasted in the excel file (delete/comment if not needed)
% clear irfsDemShk irfsSupShk;
% varplot=VarNames(:,4:end);  %selects variables only for the country-block
% irfs=zeros(NoIrfs,3,size(varplot,2));
% %PASTE the IRFs to the following spreadsheets in the corresponding xlxs file
% sheet1='IRFsDemShk86';
% sheet2='IRFsSupShk86';
% 
% 
% for i = 1:size(varplot,2)
%     if strcmpi(varplot(1,i),'GDP p. c. ')
%         varplot{i} = 'GDP';
%     end
%     if strcmpi(varplot(1,i),'10yrRate') 
%         varplot{i} = 'Long Rate';
%     end
% end
% 
% filename = strcat('../xls_files/',Country,'.xlsx');
% for i = 1:NoIrfs
%     time(i,1)=i;
% end
% VarNames1={};
% for i = 1:size(varplot,2)
%     VarNames1{1}='Quarters';
%     VarNames1(1,2+3*(i-1))=varplot(1,i);
% end
% 
% irfsDemShk=[];
% for i = 1: size(varplot,2)
%     irfs(:,:,i)=squeeze(DistDemIRFs(3+i,:,:));
%     irfsDemShk=horzcat(irfsDemShk,irfs(:,:,i));
% end
% xlswrite(filename,VarNames1,sheet1,'A1')
% xlswrite(filename,time,sheet1,'A3')
% xlswrite(filename,irfsDemShk,sheet1,'B3')
% 
% 
% irfsSupShk=[];
% for i = 1: size(varplot,2)
%     irfs(:,:,i)=squeeze(DistSuppIRFs(3+i,:,:));
%     irfsSupShk=horzcat(irfsSupShk,irfs(:,:,i));
% end
% xlswrite(filename,VarNames1,sheet2,'A1')
% xlswrite(filename,time,sheet2,'A3')
% xlswrite(filename,irfsSupShk,sheet2,'B3')