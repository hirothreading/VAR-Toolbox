%% Housekeeping
close all
clear
clc
%% Load Data
[DataMat1,XlsInfo]=xlsread('../xls_files/Dataset2.xlsx','Dataset_levels');
BeginEstim=find_variable_indices({'31/03/1975'},XlsInfo(3:end,1));
EndEstim=find_variable_indices({'31/12/2008'},XlsInfo(3:end,1));
DataMat2=DataMat1(BeginEstim:EndEstim,:);
%1 2 3 4
VarOrdeIndex=[4 7 10 11 12 13 14];
DataMat2=DataMat2(:,VarOrdeIndex);
DataMat2(:,1)=100*log(DataMat2(:,1)./100);
DataMat2(:,2)=100*log(DataMat2(:,2));
VarNames=XlsInfo(2,2:end);
VarNames=VarNames(:,VarOrdeIndex);

% DataMat3=DataMat2(2:end,:);
% DataMat3(:,1:4)=DataMat2(2:end,1:4)-DataMat2(1:end-1,1:4);

dy=cols(DataMat2);
T=rows(DataMat2);
%% SoE VAR Estimation
PriorBVARMinnesota.burn=1000;
PriorBVARMinnesota.ndraw=2000;
PriorBVARMinnesota.kskip=1;
PriorBVARMinnesota.lamda=2; %Looser (?) prior
PriorBVARMinnesota.muC=0.001;
% PriorBVARMinnesota.tau=10*PriorBVARMinnesota.lamda;
PriorBVARMinnesota.const=1;
PriorBVARMinnesota.nolags=4;
PriorBVARMinnesota.lintrend=0;
PriorBVARMinnesota.BlockIndex=[1 3;2 7];
warning off all
PostBVARMinnesota = bvar_minnesota_estim(PriorBVARMinnesota,DataMat2);
warning on all
%% VAR Identification
SignRes=zeros(dy,dy);
SignRes(1,1)=1;
SignRes(2,1)=-1;
KLbar=0;
KUbar=40;
TargetVar=1;

S=size(PostBVARMinnesota.AutoParmChain,3);
SuppIRFs=zeros(dy,12,S);
TermA=0;
for j = 1 : S
    j
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);
%     Rmat=chol(SigmaMat,'lower');
    %         Rmat=var_sign_restictions(SigmaMat,SignRes(:));%
        Rmat=uhlig(BetaMat,SigmaMat,KLbar,KUbar,TargetVar,0,0);
    if Rmat(1,1)<0;
        Rmat=-Rmat;
    end
    
    TermA=TermA+1;
    SuppIRFs(:,:,TermA)=varirf(BetaMat,Rmat(:,1),12);
end

DistSuppIRFs=quantile(SuppIRFs(:,:,1:TermA),[0.5 0.16 0.84],3);
%%
%% Plots
xxx=1:12;
figure('name','supply')
for i = 1 : dy
    subplot(3,4,i)
    plotx6(xxx,squeeze(DistSuppIRFs(i,:,:)))
    axis tight
    title(VarNames{i})
end
