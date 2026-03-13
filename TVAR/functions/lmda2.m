function out=lmda2(theta,out1,gbar,out2,sbar)

%parameter bounds
if abs(theta(1))>200 || abs(theta(2))>200
    out=1000000;
else
lamda1=theta(1);
lamda2=theta(2);


mout1=mean(out1(9:end,:),1);  %mean after two years
mout2=mean(out2(9:end,:),1);  %mean after two years
N=size(mout1,2);
% out=mean(exp(lamda1.*(mout1-repmat(gbar,1,cols(out1)))));
out=[];
for i=1:N
    GG=[mout1(:,i);mout2(:,i)];
    LL=[lamda1;lamda2]';
    MM=[gbar;sbar];
    temp=exp(LL*(GG-MM));
    out=[out;temp];
end
out=mean(out);
end
