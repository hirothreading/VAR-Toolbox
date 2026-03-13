function yhat = get_paths(N, HORZ, T, L, Y, hlast, iA, beta2, g)
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
% hlast(1)
% size(hlast)
hhat(L,:)=log(hlast(end,:));

for fi=L+1:HORZ+L
    
    xhat=[];
    for ji=1:L
        xhat=[xhat yhat(fi-ji,:)];
    end
    xhat=[xhat 1];
    
    hhat(fi,:) = hhat(fi-1,:) + ee(fi,:).*sqrt(g)';
    sigmahat   = iA * diag(exp(hhat(fi,:))) * iA';
    yhat(fi,:) = xhat*reshape(beta2(1,:),N*L+1,N) + uu(fi,:)*chol(sigmahat);
    
end
