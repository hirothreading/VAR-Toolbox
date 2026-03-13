function [ FF,MU,QQ] = comp(b,sig,n,l)



temp=reshape(b,n*l+1,n);

MU=zeros(1,n*l);
MU(1,1:n)=temp(end,:)';
FF=zeros(n*l,n*l);
FF(n+1:n*l,1:n*(l-1))=eye(n*(l-1),n*(l-1));
FF(1:n,1:n*l)=temp(1:n*l,1:n)';
QQ=zeros(n*l,n*l);
QQ(1:n,1:n)=sig;

end

