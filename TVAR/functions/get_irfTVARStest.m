function [ir,ee] = get_irfTVARStest(N, HORZIR, T, L, Y, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml)

ir=0;
ee=0;
for j=L+1:T
    y0=Y(j-L+1:j,:);
[irx,e]=get_irfTVARS4test(N, HORZIR, T, L, y0, beta1, sigma1,beta2,sigma2,tar,tarvar,tard,reps,scalex,nml);
ir=ir+irx;
ee=ee+e;
end
ir=ir/(T-L);
ee=ee/(T-L);
