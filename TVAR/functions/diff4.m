function out=diff4(x)
x1=lag0(x,4);
out=x-x1;
out=out(5:end,:);