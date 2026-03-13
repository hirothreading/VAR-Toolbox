function [out1,out2,arate]=getTgibbs(Y,X,B0,Sigma0,T0,D0,vprior,reps,burn,vscale)

out1=[];
out2=[];
T=size(X,1);
sigma2=1;
wi=ones(T,1);
omegainv=(1./sqrt(wi));

v=10;
vold=v;
naccept=0;
for i=1:reps

%step 2 Sample B conditional on sigma N(M*,V*)
xstar=X.*repmat(omegainv,1,size(X,2));
ystar=Y.*repmat(omegainv,1,size(Y,2));
XOX=xstar'*xstar;
XOY=xstar'*ystar;
V=invpd(invpd(Sigma0)+(1/sigma2)*(XOX));

M=V*(invpd(Sigma0)*B0+(1/sigma2)*XOY); 
                   
B=M+(randn(1,size(X,2))*chol(V))';


%step 2 sample sigma2
%compute residuals
resids=(ystar-xstar*B);
sigma2=ig( resids,T0,D0 );

%step 3 sample omega
resids1=Y-X*B;
temp=((resids1.^2).*(1/sigma2));
for j=1:T
    dof=v+1;
    temp1=(temp(j)+v);
    omegainv(j)=sqrt(gamm_rnd(1,1,dof/2,temp1/2));
end

%step 4 sample v using MH algorithm    
vnew=vold+randn(1,1).*sqrt(vscale);
if vnew<=0;
    accept=0;
else
   [ pnew ] = postv( vprior,vnew,omegainv.^2,T) ;
   [ pold ] = postv( vprior,vold,omegainv.^2,T) ;
   accept=exp(pnew-pold);
   
end
u=rand(1,1);
   if accept>u
       vold=vnew;
       naccept=naccept+1;
   end
v=vold;
arate=(naccept./i)*100;


if i>burn
    out1=[out1;B'];
    out2=[out2;[sigma2 v]];
end
end

end

