function ir = getirfp(N, HORZIR, T, L, Y, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml,pos)

ir=0;
parfor j=L+1:T
    y0=Y(j-L+1:j,:);
irx=getirfp4(N, HORZIR, T, L, y0, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml,pos);
ir=ir+irx;
end
ir=ir/(T-L);
