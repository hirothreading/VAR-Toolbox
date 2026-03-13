function [ omegainv,MM,VV ] = getlamda( resids1,sigma2,v )
T=rows(resids1);
temp=((resids1.^2).*(1/sigma2));
omegainv=zeros(T,1);
MM=zeros(T,1);
VV=zeros(T,1);
for j=1:T
    dof=v+1;
    temp1=(temp(j)+v);
    omegainv(j)=(gamm_rnd(1,1,dof/2,temp1/2));
    MM(j)=dof/temp1;
    VV(j)=dof;
end

end

