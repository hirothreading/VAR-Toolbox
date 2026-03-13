function [post,lik,prior]=kfilterTarnew(Y,b1,b2,hlast,iamat,L,beta0,P00,tar,tarmean,tarvariance,Ystar,ncrit)
R=0;
%%Step 2a Set up matrices for the Kalman Filter
T=rows(Y);
N=cols(Y);


beta11=beta0;
p11=P00;




x=[eye(N) zeros(N,N*(L-1))];

e1=Ystar<=tar;
if sum(e1) <ncrit ||sum(1-e1)<ncrit
    post=-inf;
    lik=-inf;
    prior=-inf;
else
% %%%%%%%%%%%Step 2b run Kalman Filter
lik=0;
for i=1:T
    H=diag(hlast(i,:));
R=iamat*H*iamat';
if e1(i)==1
    [F,mu,Q]=comp(b1,R,N,L);
else
    [F,mu,Q]=comp(b2,R,N,L);
end

    %Prediction
beta10=mu+beta11*F';
p10=F*p11*F'+Q;
%  size(beta10)
%  size(x)
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
ifeta=feta\eye(cols(feta));
K=(p10*x')*ifeta;
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

%compute likelihood
liki=-0.5*logdet(feta)+(-0.5*(eta)*ifeta*(eta'));
if isinf(liki) || ~isreal(liki) || isnan(liki)
%     lik=lik-10;
post=-inf;
    lik=-inf;
    prior=-inf;
    return
else
    lik=lik+liki;
end


end


%evaluate prior for the threshold
prior=multivariatenormal(tar,tarmean,tarvariance);

post=lik+prior;
if isinf(post) || ~isreal(post) || isnan(post)
    post=-inf;
end
end

