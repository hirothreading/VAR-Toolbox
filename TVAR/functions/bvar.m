function [bsave,ssave]=bvar(data,Lbig,REPS,BURN,Update,maxtrys,lamdaP,tauP,epsilonP)


Y=data;
L=Lbig;
N=cols(Y);
%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);

% Additional priors for VAR coefficients
muP=mean(Y)';
sigmaP=[];
deltaP=[];
e0=[];
for i=1:N
    ytemp=Y(:,i);
    xtemp=[lag0(ytemp,1) ones(rows(ytemp),1)];
    ytemp=ytemp(2:end,:);
    xtemp=xtemp(2:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/rows(ytemp);
    if abs(btemp(1))>1
        btemp(1)=1;
    end
    deltaP=[deltaP;btemp(1)];
    sigmaP=[sigmaP;stemp];
    e0=[e0 etemp];
end

[yd,xd] = create_dummiesx(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

Y0=[Y;yd];
X0=[X;xd];







T=rows(Y0);




% Display empty line (to separate samples in the screen shot)
disp(sprintf(' ')) 

  %conditional mean of the VAR coefficients
  mstar=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx=xx\eye(cols(xx));  %inv(X0'X0) to be used later in the Gibbs sampling algorithm
  sigma=eye(N); %starting value for sigma
  beta0=vec(X0\Y0);
  igibbs=1;
  jgibbs=1;
while jgibbs<REPS-BURN+1
	
    % Display progress:
    if mod(igibbs,Update)==0 
        disp(sprintf(' Replication %s of %s.', ... 
             num2str(igibbs), num2str(REPS)) );
    end
        
    
     %step 1: Sample VAR coefficients
     [ beta2,PROBLEM] = getcoef( mstar,sigma,ixx,maxtrys,N,L );
     if PROBLEM
         beta2=beta0;
     else
         beta0=beta2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma=iwpq(T,inv(scale));
    A0hat=chol(sigma);
     if igibbs>BURN && ~PROBLEM
        
      bsave(jgibbs,:)=beta2';
      ssave(jgibbs,:,:)=sigma;
jgibbs=jgibbs+1;

     end
     igibbs=igibbs+1;

end


