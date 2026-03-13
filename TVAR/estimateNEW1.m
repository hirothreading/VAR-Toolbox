clear;
beep off;

addpath('functions');
dfolder='..\data/';

REPS = 20000; %total Gibbs replications
BURN = 5000; % Burn-in
SKIP=5;  % Every skip replication is saved
Sindex=BURN+1:SKIP:REPS;
fsize=length(Sindex);
Update = 1;  %print info
MaxTrys=1000;  %
% VAR models:
VarBench.L=4;                       %lag length of the VAR
VarBench.lamdaP=0.2;                %This controls the tightness of the priors on the first lag
VarBench.tauP=10*VarBench.lamdaP;   % this controls the tightness of the priors on sum of coefficients
VarBench.epsilonP=1/1000;                % this controls tightness of the prior on the constant
VarBench.TarVAR=4;                % this specifies the column in Y that is the threshold variable
VarBench.TarD=4;                % this specifies the max delay
VarBench.TarVariance=10;                % this specifies the variance of the prior for threshold value. 
%The mean of the prior is the sample mean on line 65 (tarmean)
VarBench.TarScale=0.001;                % this specifies the variance of random walk proposal
VarBench.transform=0;  %1 for threshold variable in annual growth. 0 for no transformation
% Redefine paras (just to use shorter names):
L       = VarBench.L;
lamdaP  = VarBench.lamdaP;
tauP    = VarBench.tauP;
epsilonP= VarBench.epsilonP;
tarvar=VarBench.TarVAR;
tard=1:VarBench.TarD;
tarvariance=VarBench.TarVariance;
tarscale=VarBench.TarScale;
% Display progress:
disp(sprintf(' '));
disp(sprintf('============================================================'));
disp(sprintf('BENCHMARK MODEL: %s replications with %s burns.', num2str(REPS), num2str(BURN)));



    
% Get input data 
[databig,vnames]=xlsread('dataF.xls');
vnames=vnames(1,2:end);
data=databig;

Y=data;
N=cols(Y);
ncrit=(N*L+1);

%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

%compute threshold variable
[Ystar0,Lx]=transform(Y(:,tarvar),VarBench.transform);
Ystar=lag0(Ystar0,tard(1));


Y=Y(max([L,tard(1),Lx+tard(1)])+1:end,:);
X=X(max([L,tard(1),Lx+tard(1)])+1:end,:);
Ystar=Ystar(max([L,tard(1),Lx+tard(1)])+1:end,:);
tarmean=mean(Ystar);  %mean of the prior on the threshold is the mean value of the threshold variable
T=rows(X);
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
    sigmaP=[sigmaP;sqrt(stemp)];
    e0=[e0 etemp];
end

%dummy data to implement priors see http://ideas.repec.org/p/ecb/ecbwps/20080966.html
[yd,xd] = create_dummies(lamdaP,tauP,deltaP,epsilonP,L,muP,sigmaP,N);

save priors

 
  
  sigma1=eye(N); %starting value for sigma
  sigma2=eye(N);
  beta0=vec(X\Y);
  beta01=beta0;
  beta02=beta0;
  tar=tarmean; %initial value of the threshold 
  tarold=tar;
  
  %storage
  
   bsave1=zeros(fsize,N*(N*L+1));
   bsave2=bsave1;
   sigmaS1=zeros(fsize,N,N);
   sigmaS2=zeros(fsize,N,N);  
   tsave=zeros(fsize,2);
   regime=zeros(fsize,T);
  
  naccept=0;
  igibbs=1;
  jgibbs=1;
while jgibbs<=fsize
	
    
    
    %step 1: Seperate into two regimes
    e1=Ystar<=tar;
    e2=Ystar>tar;
    
    Y1=Y(e1,:);
    X1=X(e1,:);
    
    Y2=Y(e2,:);
    X2=X(e2,:);
    
    %step 2 Sample Coefficients and variance regime 1
    
    Y0=[Y1;yd];
    X0=[X1;xd];
  %conditional mean of the VAR coefficients
  mstar1=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx1=xx\eye(cols(xx));
   [ beta1,PROBLEM1] = getcoef( mstar1,sigma1,ixx1,MaxTrys,N,L );
     if PROBLEM1
         beta1=beta01;
     else
         beta01=beta1;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta1,N*L+1,N);
    scale=e'*e;
    sigma1=iwpQ(rows(Y0),inv(scale));  
    
    
    %step 3 Sample Coefficients and variance in regime 2
    
     Y0=[Y2;yd];
    X0=[X2;xd];
  %conditional mean of the VAR coefficients
  mstar2=vec(X0\Y0);  %ols on the appended data
  xx=X0'*X0;
  ixx2=xx\eye(cols(xx));
   [ beta2,PROBLEM2] = getcoef( mstar2,sigma2,ixx2,MaxTrys,N,L );
     if PROBLEM2
         beta2=beta02;
     else
         beta02=beta2;
     end
     
     %draw covariance
     e=Y0-X0*reshape(beta2,N*L+1,N);
    scale=e'*e;
    sigma2=iwpQ(rows(Y0),inv(scale)); 
    
    
    
    %step 4 Sample Threshold via a Random Walk Metropolis Step
    
    tarnew=tarold+randn(1,1)*sqrt(tarscale);
     %compute conditional posterior at the old and new draw

      postnew=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarnew,tarmean,tarvariance,Ystar,ncrit);
      postold=getvarpost(Y,X,beta1,beta2,sigma1,sigma2,L,tarold,tarmean,tarvariance,Ystar,ncrit);

     accept=exp(postnew-postold);
     u=rand(1,1);
     if u<accept
         tarold=tarnew;
         naccept=naccept+1;
     end
     tar=tarold;
     arate=naccept/igibbs;
     
     %step 5 sample delay parameter
     prob=[];
     for jj=1:length(tard)
         [yy,xx,yys]=preparexx( data,L,tard(jj),tarvar,VarBench );
         
         postx=getvarpost(yy,xx,beta1,beta2,sigma1,sigma2,L,tar,tarmean,tarvariance,yys,ncrit);
         prob=[prob postx];
     end
     prob=safeexp(prob);
     prob=prob/sum(prob);
     
     %draw delay
     index=discretesample(prob,1);
     %update data
     [Y,X,Ystar]=preparexx( data,L,tard(index),tarvar,VarBench );
     
     
     
    % Display progress:
    if mod(igibbs,Update)==0 
        disp(sprintf(' Sample %s    Replication %s of %s acceptance %s. delay is %s.', ... 
            num2str(1), num2str(igibbs), num2str(REPS),num2str(arate),num2str(tard(index)) ));
    end 
     
     
     
%      if igibbs>100 && igibbs<5000
%          if arate<0.35
%              tarscale=tarscale*0.99;
%          elseif arate>0.55
%              tarscale=tarscale*1.01;
%          end
%      end
    
    

     if igibbs>BURN 
         if sum(Sindex==igibbs)>0
        





   bsave1(jgibbs,:)=beta1';
   bsave2(jgibbs,:)=beta2';
   sigmaS1(jgibbs,:,:)=sigma1;
   sigmaS2(jgibbs,:,:)=sigma2;   
   tsave(jgibbs,:)=[tar tard(index)];
    jgibbs=jgibbs+1;
         end
     end
     igibbs=igibbs+1;
end




% Display progress:
disp(sprintf(' '));
disp(sprintf('BENCHMARK MODEL - done'));

save('results','bsave1','bsave2','sigmaS1','sigmaS2','tsave','regime')

%plot of threshold variable

figure(1)
tmp=mean(tsave);
     [Y,X,Ystar]=preparexx( data,L,round(tmp(2)),tarvar,VarBench );
    e1=Ystar<=tar;

  TT=1:rows(e1);
  names={'S_{t}<=Y^{*}';'Threshold Variable'};
   plotregx( TT,e1,Ystar,names )
 
  title('Threshold')
