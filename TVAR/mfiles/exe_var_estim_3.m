%% Housekeeping
close all
clear
clc
%% Load Data
[DataMat1,XlsInfo]=xlsread('../xls_files/Dataset3.xlsx','Dataset_levels');
BeginEstim=find_variable_indices({'31/03/1975'},XlsInfo(3:end,1));
EndEstim=find_variable_indices({'31/12/2015'},XlsInfo(3:end,1));
DataMat2=DataMat1(BeginEstim:EndEstim,:);
VarNames=XlsInfo(2,2:end);
VarNames=VarNames(:,[4 7 9 10 11 12 13 14]);
DataMat=[100*log(DataMat2(:,4)./100) 100*log(DataMat2(:,8)./DataMat2(:,2))...
    DataMat2(:,9:end)];
dy=cols(DataMat);
T=rows(DataMat);
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
PriorBVARMinnesota.BlockIndex=[1 3;2 8];
warning off all
PostBVARMinnesota = bvar_minnesota_estim(PriorBVARMinnesota,DataMat);
warning on all
%% VAR Identification
S=size(PostBVARMinnesota.AutoParmChain,3);

pp = {'N:\Cross Divisional Work\DSGE Applications Projects\EnergyProjectAhmed\m_files'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\Bayesianindirict',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\DistStats',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\TestStat',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\MAModels',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\HACC',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\dsge_estim',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\DSGEANALYSIS'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\func',...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\sims_Optimization'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\solution'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\solution\SP_SOLVE'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\VAR'...
    'N:\Cross Divisional Work\DSGE Applications Projects\ToolKit\VAR\BVAR'};

jobManagerProfile=parallel.importProfile('../../MATLAB_HEAD.settings');
jobManager = parcluster(jobManagerProfile);
job = createJob(jobManager);
ncpu = 32;
% set(job,'NumWorkersRange',ncpu);
disp(['Running on matlab cluster using ',num2str(ncpu),' workers']);

SignRes=zeros(dy,dy);


for j = 1 : S
    BetaMat=PostBVARMinnesota.AutoParmChain(:,PriorBVARMinnesota.const+PriorBVARMinnesota.lintrend+1:end,j);
    SigmaMat=PostBVARMinnesota.CovParmChain(:,:,j);

    createTask(job,@median_var_sign_penalty_function,3,...
        {SigmaMat,BetaMat,12,4});
end

disp('Job submitted!');
alltasks = get(job,'Tasks');
taskErrors = get(alltasks,'ErrorMessage');
set(job,'AdditionalPaths',pp,'AttachedFiles',{})
set(alltasks,'Timeout',60,'MaximumRetries',0)
submit(job);                                 % submit tasks
wait(job,'finished');                           % wait until finished
[NewOutputs,IndexSucces]=collect_parallel_output(alltasks,3);
get(job,'StartTime')
get(job,'FinishTime')
delete(job);

SuppIRFs=zeros(dy,12,rows(NewOutputs));
TermA=0;

for i = 1 : rows(NewOutputs);
    SuppIRFs(:,:,i)=NewOutputs{i,2};
end

DistSuppIRFs=quantile(SuppIRFs,[0.5 0.16 0.84],3);
%% Plots
xxx=1:12;
figure('name','supply')
for i = 1 : dy
    subplot(3,4,i)
    plotx6(xxx,squeeze(DistSuppIRFs(i,:,:)))
    axis tight
    title(VarNames{i})
end

