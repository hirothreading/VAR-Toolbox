function [ ir1,irp1,yy,yy1,p,p1] = getimpulse(N, HORZ, L, Y0, beta1, ...
    sigma1,beta2,sigma2,tar,tvar,delay,reps,Lx,A01,A02,Y0m,scale,pos,transform)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


csigma1=chol(sigma1);
csigma2=chol(sigma2);
b1x=reshape(beta1,N*L+1,N);
b2x=reshape(beta2,N*L+1,N);
LL=max([L Lx+delay]);


yy1=0;
yy=0;
p1=0;
p=0;
for ii=1:reps
%compute forecast
yhat=zeros(HORZ+LL,N);
yhat(1:LL,:)=Y0;
yhat1=yhat;


yhatm=zeros(HORZ+LL,N);
yhatm(1:LL,:)=Y0m;
yhat1m=yhatm;


ystar=zeros(HORZ+LL,1);
ystar1=ystar;



pstar=zeros(HORZ+LL,1);
pstar1=ystar;





for fi=LL+1:HORZ+LL
    
    
    xhat=ones(1,N*L+1);
    xhat1=ones(1,N*L+1);
    
    
    jix=1;
    for ji=1:L
        xhat(:,jix:jix+N-1)=yhat(fi-ji,:);
         xhat1(:,jix:jix+N-1)=yhat1(fi-ji,:);
        
       
         jix=jix+N;
    end
  
    if transform==0
    yhatm(fi,:)=yhat(fi,:);
    yhat1m(fi,:)=yhat1(fi,:);
    elseif transform==1
    
       % [fi yhat1(fi,:) yhat1(fi-Lx,:)]
    yhatm(fi,:)=yhat(fi-delay,:)-yhat(fi-Lx-delay,:);
    yhat1m(fi,:)=yhat1(fi-delay,:)-yhat1(fi-Lx-delay,:);
    end
        

    ystar(fi,:)=(yhatm(fi,tvar));
    ystar1(fi,:)=yhat1m(fi,tvar);
   

    
    
    %simulate data
    if fi==LL+1
    uu1=zeros(1,N);uu2=zeros(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*b1x + uu1*csigma1)*e1)+...
        ((xhat*b2x + uu2*csigma2)*e2);
    pstar(fi,:)=e1;
    %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(pos)=scale;uu2(pos)=scale;
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    
    yhat1(fi,:) = ((xhat1*b1x + uu1*A01)*e1)+...
        ((xhat1*b2x + uu2*A02)*e2);
     pstar1(fi,:)=e1;
     
    else
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*b1x + uu1*csigma1)*e1)+...
        ((xhat*b2x + uu2*csigma2)*e2);
        pstar(fi,:)=e1;
    
    %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    yhat1(fi,:) = ((xhat1*b1x + uu1*csigma1)*e1)+...
        ((xhat1*b2x + uu2*csigma2)*e2);
        pstar1(fi,:)=e1;
     
   
    end
   
end

yy=yy+yhat;
  yy1=yy1+yhat1;

 
  
  
  p=p+pstar;
   
  p1=p1+pstar1;
 
 
end

yy=yy/reps;
yy1=yy1/reps;



p=p/reps;
p1=p1/reps;








ir1=yy1-yy;

ir1=ir1(LL+1:end,:);



irp1=p1-p;
irp1=irp1(LL+1:end,:);

