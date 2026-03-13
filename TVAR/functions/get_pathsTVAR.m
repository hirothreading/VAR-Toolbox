function [yhat ystar] = ...
    get_pathsTVAR(N, HORZ, T, L, Y, beta1, sigma1,beta2,sigma2,tar,tvar,delay)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


% Draw N(0,1) innovations for variance and mean equation:
csigma1=chol(sigma1);
csigma2=chol(sigma2);

uu1 = randn(HORZ+L,N);
uu2 = randn(HORZ+L,N);
% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.

%compute forecast
yhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y(T-L+1:T,:);
ystar=zeros(HORZ+L,1);
for fi=L+1:HORZ+L
    
    xhat=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
    end
    xhat=[xhat 1];
   
    ystar(fi,:)=yhat(fi-delay,tvar);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = ((xhat*reshape(beta1,N*L+1,N) + uu1(fi,:)*csigma1)*e1)+...
        ((xhat*reshape(beta2,N*L+1,N) + uu2(fi,:)*csigma2)*e2);
    
end
