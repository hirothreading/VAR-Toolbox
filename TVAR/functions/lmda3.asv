function out=lmda3(theta,out1,gbar,out2,sbar)

%parameter bounds
check=abs(theta)>200;
if sum(check)>0
    out=1000000;
else
LL=theta';
mout1=out1;
mout2=out2;
N=size(mout1,2);
% out=mean(exp(lamda1.*(mout1-repmat(gbar,1,cols(out1)))));
out=[];
for i=1:N
    GG=[mout1(:,i);mout2(:,i)];
    MM=[gbar;sbar];
    temp=exp(LL*(GG-MM));
    out=[out;temp];
end
out=mean(out);
end
