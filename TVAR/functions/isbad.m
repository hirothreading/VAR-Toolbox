function [ check ] = isbad( x )
check1=isnan(x);
check2=1-isreal(x);
check3=isinf(x);

check=(check1+check2+check3)>0;
end

