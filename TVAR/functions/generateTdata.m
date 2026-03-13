function y=generateTdata(beta,sigmasquared,nobs,DF,X)
mu=X*beta;
%generate y from the student tdensity
y=trnd(DF,nobs,1).*sqrt(sigmasquared)+mu;
