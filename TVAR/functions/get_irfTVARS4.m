function [ ir4,yy4,yy] = get_irfTVARS4(N, HORZ, T, L, Y0, beta1, sigma1,beta2,sigma2,tar,tvar,delay,reps,scale,nml)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


% Draw N(0,1) innovations for variance and mean equation:
csigma1=chol(sigma1);
csigma2=chol(sigma2);
d1=diag(csigma1);
d2=diag(csigma2);
csx1=csigma1./repmat(d1,1,N);
csx2=csigma2./repmat(d2,1,N);

% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.
yy1=0;
yy2=0;
yy3=0;yy4=0;yy=0;
for ii=1:reps
%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y0;

yhat4=yhat;
ystar=zeros(HORZ+L,1);

ystar4=ystar;

for fi=L+1:HORZ+L
    
    
    xhat=[];
    xhat1=[];
    xhat2=[];
    xhat3=[];
    xhat4=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
        xhat4=[xhat4 yhat4(fi-ji,:)];
    end
    xhat=[xhat 1];xhat4=[xhat4 1];
  
    ystar(fi,:)=yhat(fi-delay,tvar);
    ystar4(fi,:)=yhat4(fi-delay,tvar);
    
    
    %simulate data
    if fi==L+1
    uu1=zeros(1,N);uu2=zeros(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    
     %
     
    uu1=zeros(1,N);uu2=zeros(1,N);uu1(4)=scale;uu2(4)=scale;
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    if nml==0
    yhat4(fi,:) = ((xhat4*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat4*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    else
        yhat4(fi,:) = ((xhat4*reshape(beta1,N*L+1,N) + uu1*csx1)*e1)+...
        ((xhat4*reshape(beta2,N*L+1,N) + uu2*csx2)*e2);
    end
    else
        uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    
     %
    uu1=randn(1,N);uu2=randn(1,N);
    e1=ystar4(fi,:)<=tar;
    e2=ystar4(fi,:)>tar;
    yhat4(fi,:) = ((xhat4*reshape(beta1,N*L+1,N) + uu1*csigma1)*e1)+...
        ((xhat4*reshape(beta2,N*L+1,N) + uu2*csigma2)*e2);
    end
   
end

yy=yy+yhat;
   
   yy4=yy4+yhat4;

end

yy=yy/reps;

yy4=yy4/reps;


ir4=yy4-yy;


 ir4=ir4(L+1:end,:);
