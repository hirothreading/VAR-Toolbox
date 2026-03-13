function out=lmda(theta,out1,gbar)

%parameter bounds
if abs(theta)>200
    out=1000000;
else
lamda1=theta;

mout1=mean(out1,1);
out=mean(exp(lamda1.*(mout1-repmat(gbar,1,cols(out1)))));
end
