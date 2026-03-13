function Y=gampdf2(mu,v,h)
X=h;
A=v/2;
B=(2*mu)/v;

Y = log(gampdf(X,A,B));