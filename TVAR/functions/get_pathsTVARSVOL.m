function yhat = ...
    get_pathsTVARSVOL(N, HORZ, T, L, Y, hlast, iA, beta1,beta2, g,tar,tvar,delay)
% -------------------------------------------------------------------------
% get_paths:
% generates a matrix of simulated paths for Y for given parameters and
% general set-up (lags, horizon, etc.)
% -------------------------------------------------------------------------


% Draw N(0,1) innovations for variance and mean equation:
ee = randn(HORZ+L,N);
uu = randn(HORZ+L,N);
% Note we only need HORZ*N innovations, but adding an extra L draws makes 
% the indexing in the loop below a bit cleaner.

%compute forecast
yhat=zeros(HORZ+L,N);
hhat=zeros(HORZ+L,N);
yhat(1:L,:)=Y(T-L+1:T,:);
hhat(L,:)=log(hlast(end,:));
ystar=zeros(HORZ+L,1);

for fi=L+1:HORZ+L
    
    xhat=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
    end
    xhat=[xhat 1];
    
    hhat(fi,:) = hhat(fi-1,:) + ee(fi,:).*sqrt(g);
    sigmahat   = iA * diag(exp(hhat(fi,:))) * iA';
    ystar(fi,:)=yhat(fi-delay,tvar);
    e1=ystar(fi,:)<=tar;
    e2=ystar(fi,:)>tar;
    yhat(fi,:) = (xhat*reshape(beta1(1,:),N*L+1,N)).*e1 +...
        (xhat*reshape(beta2(1,:),N*L+1,N)).*e2 +...
        uu(fi,:)*cholx(sigmahat);
    
end
