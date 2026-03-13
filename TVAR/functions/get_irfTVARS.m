function ir = get_irfTVARS(N, HORZIR, T, L, Y, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml)

ir=0;
parfor j=L+1:T
    y0=Y(j-L+1:j,:);
irx=get_irfTVARS4(N, HORZIR, T, L, y0, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml);
ir=ir+irx;
end
ir=ir/(T-L);
