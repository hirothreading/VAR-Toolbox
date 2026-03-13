%%% 1. This file plots the IRFs for different countries in response to
% oil demand/supply shocks
%%% 2. Inputs requried: IRFs for each country must be pasted into the excel
% files before running this code. These IRFs are created in
% exe_var_estim_signs.m (see: IRFsDemShk and IRFsSupShk)

%% Housekeeping
%close all
clear
clc

%% Specify Countries and Variables
Country={'US','Germany','France','Spain','Italy'};%,'Australia','Canada','Japan','Switzerland','Italy'};
VarNames={'GDP','CPI','PPI','Policy Rate','Long Rate'};
NoIrfs=20;

%% Plot Figure
figure%('name',Country{i})
for i = 1 : size(Country,2)
    [DataMat1,XlsInfo]=xlsread(strcat('../xls_files/',Country{i},'.xlsx'),'IRFsSupShk86');
    
    for j = 1 : size(VarNames,2)
        subplot(size(Country,2),size(VarNames,2),(i-1)*size(VarNames,2)+j)
                
        VarNo=find_variable_indices({VarNames(1,j)},XlsInfo(1,1:end))        
        if VarNo == 0
            if i == 1
                axis tight
                title(VarNames{j})
            end   
            tx = text(0.20,0.5,'Insufficient Data');
            set(tx,'fontweight','bold');
            continue
        end
        
        plot(DataMat1(1:NoIrfs,1),DataMat1(1:NoIrfs,VarNo),'k',DataMat1(1:NoIrfs,1),DataMat1(1:NoIrfs,VarNo+1),'-.r',DataMat1(1:NoIrfs,1),DataMat1(1:NoIrfs,VarNo+2),'-.r','linewidth',1.5)
        hline = refline([0 0]);
        hline.Color = 'b';
        axis tight
        if j == 1
            ylabel(Country{i},'fontweight','bold','fontsize',13)
        end
        if i == 1
            title(VarNames{j},'fontsize',13)
        end
        if i == size(Country,2) && j == 1 %size(VarNames,2)
            legend('Median','Conf. Band','location','southwest')
            legend('boxoff')
        end
   end
end