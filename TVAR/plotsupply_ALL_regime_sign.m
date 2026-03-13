clear
addpath('functions')
shock=1;
PP=[50 16 84];


load('priors')


load irf1sd_r1_sign
irf1_r1=irf;

%
load irf1sd_r2_sign
irf1_r2=irf;




HH=(0:39)';
 figure(1)


for variable=1:N



out1_r1=prctile(irf1_r1(:,:,variable),PP);
out1_r2=prctile(irf1_r2(:,:,variable),PP);






 subplot(2,2,variable)
 [~,hx]=plotx2(HH,out1_r1');
 hold on
 h=plot(HH,out1_r2',':','color',[0.5 0.1 0.4]);
 set(h(1),'LineWidth',5.5);
 set(h(2),'LineWidth',1.5);
 set(h(3),'LineWidth',1.5);

 
xlim([0 30])

 title(vnames{variable},'interpreter','latex','FontSize',16)
 
 xlabel('Horizon')
 
 if variable==N
      legend([hx  h(1)],'Regime 1','Regime 2')
 end
end

