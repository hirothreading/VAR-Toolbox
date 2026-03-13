function [ ir1,ir2,ir3,ir4,ir5,ir6,ir7,irp1,irp2,irp3,irp4,irp5,irp6,irp7] = getimpulse7(N, HORZ, L, Y0, beta1, ...
    sigma1,beta2,sigma2,tar,tvar,delay,reps,Lx,A01,A02,Y0m,scale)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


csigma1=chol(sigma1);
csigma2=chol(sigma2);
b1x=reshape(beta1,N*L+1,N);
b2x=reshape(beta2,N*L+1,N);


yy1=0;
yy2=0;
yy3=0;yy4=0;yy5=0;yy6=0;yy7=0;yy=0;
p1=0;
p2=0;
p3=0;p4=0;p5=0;p6=0;p7=0;p=0;
for ii=1:reps
%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y0;
yhat1=yhat;
yhat2=yhat;
yhat3=yhat;
yhat4=yhat;
yhat5=yhat;
yhat6=yhat;
yhat7=yhat;
yhatm=zeros(HORZ+L,N);
yhatm(1:L,:)=Y0m;
yhat1m=yhatm;
yhat2m=yhatm;
yhat3m=yhatm;
yhat4m=yhatm;
yhat5m=yhatm;
yhat6m=yhatm;
yhat7m=yhatm;
ystar=zeros(HORZ+L,1);
ystar1=ystar;
ystar2=ystar;
ystar3=ystar;
ystar4=ystar;
ystar5=ystar;
ystar6=ystar;
ystar7=ystar;

pstar=zeros(HORZ+L,1);
pstar1=ystar;
pstar2=ystar;
pstar3=ystar;
pstar4=ystar;
pstar5=ystar;
pstar6=ystar;
pstar7=ystar;




for fi=L+1:HORZ+L
    
    
    xhat=ones(1,N*L+1);
    xhat1=ones(1,N*L+1);
    xhat2=ones(1,N*L+1);
    xhat3=ones(1,N*L+1);
    xhat4=ones(1,N*L+1);
    xhat5=ones(1,N*L+1);
    xhat6=ones(1,N*L+1);
    xhat7=ones(1,N*L+1);
    jix=1;
    for ji=1:L
        xhat(:,jix:jix+N-1)=yhat(fi-ji,:);
         xhat1(:,jix:jix+N-1)=yhat1(fi-ji,:);
         xhat2(:,jix:jix+N-1)=yhat2(fi-ji,:);
       xhat3(:,jix:jix+N-1)=yhat3(fi-ji,:);
        xhat4(:,jix:jix+N-1)=yhat4(fi-ji,:);
        xhat5(:,jix:jix+N-1)=yhat5(fi-ji,:);
         xhat6(:,jix:jix+N-1)=yhat6(fi-ji,:);
         xhat7(:,jix:jix+N-1)=yhat7(fi-ji,:);
         jix=jix+N;
    end
  
    
    yhatm(fi,:)=yhat(fi,:);
    yhat1m(fi,:)=yhat1(fi,:);
     yhat2m(fi,:)=yhat2(fi,:);
      yhat3m(fi,:)=yhat3(fi,:);
       yhat4m(fi,:)=yhat4(fi,:);
        yhat5m(fi,:)=yhat5(fi,:);
         yhat6m(fi,:)=yhat6(fi,:);
          yhat7m(fi,:)=yhat7(fi,:);
    for jj=1:Lx-1
        yhatm(fi,:)= yhatm(fi,:)+yhat(fi-jj,:);
        yhat1m(fi,:)= yhat1m(fi,:)+yhat1(fi-jj,:);
        yhat2m(fi,:)= yhat2m(fi,:)+yhat2(fi-jj,:);
        yhat3m(fi,:)= yhat3m(fi,:)+yhat3(fi-jj,:);
        yhat4m(fi,:)= yhat4m(fi,:)+yhat4(fi-jj,:);
        yhat5m(fi,:)= yhat5m(fi,:)+yhat5(fi-jj,:);
        yhat6m(fi,:)= yhat6m(fi,:)+yhat6(fi-jj,:);
        yhat7m(fi,:)= yhat7m(fi,:)+yhat7(fi-jj,:);
    end
    ystar(fi,:)=yhatm(fi-delay,tvar);
    ystar1(fi,:)=yhat1m(fi-delay,tvar);
    ystar2(fi,:)=yhat2m(fi-delay,tvar);
    ystar3(fi,:)=yhat3m(fi-delay,tvar);
    ystar4(fi,:)=yhat4m(fi-delay,tvar);
    ystar5(fi,:)=yhat5m(fi-delay,tvar);
    ystar6(fi,:)=yhat6m(fi-delay,tvar);
    ystar7(fi,:)=yhat7m(fi-delay,tvar);
    
    
    %simulate data
    if fi==L+1
    uu1=zeros(1,N);uu2=zeros(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*b1x + uu1*csigma1)*e1)+...
        ((xhat*b2x + uu2*csigma2)*e2);
    pstar(fi,:)=e1;
    %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(1)=scale;uu2(1)=scale;
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    
    yhat1(fi,:) = ((xhat1*b1x + uu1*A01)*e1)+...
        ((xhat1*b2x + uu2*A02)*e2);
     pstar1(fi,:)=e1;
       %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(2)=scale;uu2(2)=scale;
    e1=ystar2(fi,:)<=tar;
    e2=ystar2(fi,:)>tar;
    
    yhat2(fi,:) = ((xhat2*b1x + uu1*A01)*e1)+...
        ((xhat2*b2x + uu2*A02)*e2);
        pstar2(fi,:)=e1;
       %
    
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(3)=scale;uu2(3)=scale;
    e1=ystar3(fi,:)<=tar;
    e2=ystar3(fi,:)>tar;
    
    yhat3(fi,:) = ((xhat3*b1x + uu1*A01)*e1)+...
        ((xhat3*b2x + uu2*A02)*e2);
        pstar3(fi,:)=e1;
     %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(5)=scale;uu2(5)=scale;
    e1=ystar5(fi,:)<=tar;
    e2=ystar5(fi,:)>tar;
    yhat5(fi,:) = ((xhat5*b1x + uu1*A01)*e1)+...
        ((xhat5*b2x + uu2*A02)*e2);
        pstar5(fi,:)=e1;
     %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(6)=scale;uu2(6)=scale;
    e1=ystar6(fi,:)<=tar;
    e2=ystar6(fi,:)>tar;
    
    yhat6(fi,:) = ((xhat6*b1x + uu1*A01)*e1)+...
        ((xhat6*b2x + uu2*A02)*e2);
        pstar6(fi,:)=e1;
        %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(7)=scale;uu2(7)=scale;
    e1=ystar7(fi,:)<=tar;
    e2=ystar7(fi,:)>tar;
    
    yhat7(fi,:) = ((xhat7*b1x + uu1*A01)*e1)+...
        ((xhat7*b2x + uu2*A02)*e2);
        pstar7(fi,:)=e1;
    
    
    
     %
     
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(4)=scale;uu2(4)=scale;
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    
    yhat4(fi,:) = ((xhat4*b1x + uu1*A01)*e1)+...
        ((xhat4*b2x + uu2*A02)*e2);
        pstar4(fi,:)=e1;
    else
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*b1x + uu1*csigma1)*e1)+...
        ((xhat*b2x + uu2*csigma2)*e2);
        pstar(fi,:)=e1;
     %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    yhat4(fi,:) = ((xhat4*b1x + uu1*csigma1)*e1)+...
        ((xhat4*b2x + uu2*csigma2)*e2);
        pstar4(fi,:)=e1;
    %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    yhat1(fi,:) = ((xhat1*b1x + uu1*csigma1)*e1)+...
        ((xhat1*b2x + uu2*csigma2)*e2);
        pstar1(fi,:)=e1;
      %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar2(fi,:)<=tar;
    e2=ystar2(fi,:)>tar;
    
    yhat2(fi,:) = ((xhat2*b1x + uu1*csigma1)*e1)+...
        ((xhat2*b2x + uu2*csigma2)*e2);
        pstar2(fi,:)=e1;
       %
   uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar3(fi,:)<=tar;
    e2=ystar3(fi,:)>tar;
    yhat3(fi,:) = ((xhat3*b1x + uu1*csigma1)*e1)+...
        ((xhat3*b2x + uu2*csigma2)*e2);
        pstar3(fi,:)=e1;
     %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar5(fi,:)<=tar;
    e2=ystar5(fi,:)>tar;
    yhat5(fi,:) = ((xhat5*b1x + uu1*csigma1)*e1)+...
        ((xhat5*b2x + uu2*csigma2)*e2);
        pstar5(fi,:)=e1;
     %
   uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar6(fi,:)<=tar;
    e2=ystar6(fi,:)>tar;
    
    yhat6(fi,:) = ((xhat6*b1x + uu1*csigma1)*e1)+...
        ((xhat6*b2x + uu2*csigma2)*e2);
        pstar6(fi,:)=e1;
        %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar7(fi,:)<=tar;
    e2=ystar7(fi,:)>tar;
    
    yhat7(fi,:) = ((xhat7*b1x + uu1*csigma1)*e1)+...
        ((xhat7*b2x + uu2*csigma2)*e2);
        pstar7(fi,:)=e1;
    
    
    
    
    
    
    end
   
end

yy=yy+yhat;
   yy4=yy4+yhat4;
  yy1=yy1+yhat1;
  yy2=yy2+yhat2;
  yy3=yy3+yhat3;
  yy5=yy5+yhat5;
  yy6=yy6+yhat6;
  yy7=yy7+yhat7;
  
  
  p=p+pstar;
   p4=p4+pstar4;
  p1=p1+pstar1;
  p2=p2+pstar2;
  p3=p3+pstar3;
  p5=p5+pstar5;
  p6=p6+pstar6;
  p7=p7+pstar7;
end

yy=yy/reps;
yy4=yy4/reps;
yy1=yy1/reps;
yy2=yy2/reps;
yy3=yy3/reps;
yy5=yy5/reps;
yy6=yy6/reps;
yy7=yy7/reps;

p=p/reps;
p1=p1/reps;p2=p2/reps;p3=p3/reps;p4=p4/reps;
p5=p5/reps;p6=p6/reps;p7=p7/reps;








ir1=yy1-yy;
ir2=yy2-yy;
ir3=yy3-yy;
ir5=yy5-yy;
ir6=yy6-yy;
ir4=yy4-yy;
ir7=yy7-yy;
ir4=ir4(L+1:end,:);
ir1=ir1(L+1:end,:);
ir2=ir2(L+1:end,:);
ir3=ir3(L+1:end,:);
ir5=ir5(L+1:end,:);
ir6=ir6(L+1:end,:);
ir7=ir7(L+1:end,:);



irp1=p1-p;
irp2=p2-p;
irp3=p3-p;
irp5=p5-p;
irp6=p6-p;
irp4=p4-p;
irp7=p7-p;
irp4=irp4(L+1:end,:);
irp1=irp1(L+1:end,:);
irp2=irp2(L+1:end,:);
irp3=irp3(L+1:end,:);
irp5=irp5(L+1:end,:);
irp6=irp6(L+1:end,:);
irp7=irp7(L+1:end,:);