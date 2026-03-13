function [hh,hh1]=plotx2(t,y)
set(gcf,'DefaultAxesColorOrder',[0.8 0.1 0.1;1 0 0;1 0 0;0 0 1]);
cu=y(:,2);
cl=y(:,3);

h=t;
if cols(h)>1
h=h';
end
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[1 0.9 0.9]);
set(hh,'facecolor',[1 0.9 0.9]);


hold on
hh1=plot(h,y(:,1),'r','LineWidth',5.5);
axis tight
 hold on;
zz=zeros(size(y,1),1);
 plot(h,zz,'b-');
