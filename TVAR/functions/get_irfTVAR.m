function [ir1, ir2, ir3, ir4] = get_irfTVAR(N, HORZ, T, L, Y, beta1, sigma1,beta2,sigma2,tar,tvar,delay,reps)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


% Draw N(0,1) innovations for variance and mean equation:
csigma1=chol(sigma1);
csigma2=chol(sigma2);

% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.
yy1=0;
yy2=0;
yy3=0;yy4=0;yy=0;
for ii=1:reps
%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y(T-L+1:T,:);
yhat1=yhat;
yhat2=yhat;
yhat3=yhat;
yhat4=yhat;
ystar=zeros(HORZ+L,1);
ystar1=ystar;
ystar2=ystar;
ystar3=ystar;
ystar4=ystar;

for fi=L+1:HORZ+L
    
    
    xhat=[];
    xhat1=[];
    xhat2=[];
    xhat3=[];
    xhat4=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
        xhat1=[xhat1 yhat1(fi-ji,:)];
        xhat2=[xhat2 yhat2(fi-ji,:)];
        xhat3=[xhat3 yhat3(fi-ji,:)];
        xhat4=[xhat4 yhat4(fi-ji,:)];
    end
    xhat=[xhat 1];xhat1=[xhat1 1];xhat2=[xhat2 1];xhat3=[xhat3 1];xhat4=[xhat4 1];
  
    ystar(fi,:)=yhat(fi-delay,tvar);
    ystar1(fi,:)=yhat1(fi-delay,tvar);
    ystar2(fi,:)=yhat2(fi-delay,tvar);
    ystar3(fi,:)=yhat3(fi-delay,tvar);
    ystar4(fi,:)=yhat4(fi-delay,tvar);
    
    
    %simulate data
    if fi==L+1
    uu1=zeros(1,N);uu2=zeros(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(1)=1;uu2(1)=1;
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    yhat1(fi,:) = ((xhat1*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat1*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    
    %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(2)=1;uu2(2)=1;
    e1=ystar2(fi,:)<=tar;
    e2=ystar2(fi,:)>tar;
    yhat2(fi,:) = ((xhat2*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat2*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
%     yhat2-yhat
     %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(3)=1;uu2(3)=1;
    e1=ystar3(fi,:)<=tar;
    e2=ystar3(fi,:)>tar;
    yhat3(fi,:) = ((xhat3*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat3*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
     %
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(4)=1;uu2(4)=1;
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    yhat4(fi,:) = ((xhat4*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat4*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    else
        uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar1(fi,:)<=tar;
    e2=ystar1(fi,:)>tar;
    yhat1(fi,:) = ((xhat1*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat1*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    %
   uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar2(fi,:)<=tar;
    e2=ystar2(fi,:)>tar;
    yhat2(fi,:) = ((xhat2*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat2*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
     %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar3(fi,:)<=tar;
    e2=ystar3(fi,:)>tar;
    yhat3(fi,:) = ((xhat3*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat3*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
     %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    yhat4(fi,:) = ((xhat4*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat4*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    end
   yy=yy+yhat;
   yy1=yy1+yhat1;
   yy2=yy2+yhat2;
   yy3=yy3+yhat3;
   yy4=yy4+yhat4;
end
yy=yy/reps;
yy1=yy1/reps;
yy2=yy2/reps;
yy3=yy3/reps;
yy4=yy4/reps;

ir1=yy1-yy;
ir2=yy2-yy;
ir3=yy3-yy;
ir4=yy4-yy;

ir1=ir1(L+1:end,:);
 ir2=ir2(L+1:end,:);   
 ir3=ir3(L+1:end,:);
 ir4=ir4(L+1:end,:);
end
