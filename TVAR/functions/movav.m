function out=movav(x,L)
out=x;
for j=1:L
    out=out+lag0(x,j);
end
out=out./(L+1);