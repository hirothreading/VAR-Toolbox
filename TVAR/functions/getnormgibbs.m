function [out1,out2]=getnormgibbs(Y,X,B0,Sigma0,T0,D0,reps,burn)

out1=[];
out2=[];
T=size(X,1);
sigma2=1;
for i=1:reps

%step 2 Sample B conditional on sigma N(M*,V*)
M=inv(inv(Sigma0)+(1/sigma2)*(X'*X))*(inv(Sigma0)*B0+(1/sigma2)*X'*Y); 
V=inv(inv(Sigma0)+(1/sigma2)*(X'*X));
                   
B=M+(randn(1,size(X,2))*chol(V))';


%step 3 sample sigma2 conditional on B from IG(T1,D1);
%compute residuals
resids=Y-X*B;
%compute posterior df and scale matrix
T1=T0+T;
D1=D0+resids'*resids;
%draw from IG
z0=randn(T1,1);
z0z0=z0'*z0;
sigma2=D1/z0z0;

if i>burn
    out1=[out1;B'];
    out2=[out2;sigma2];
end
end

end

