function out=movsum(x,L)
out=x;
for j=1:L
    out=out+lag0(x,j);
end
