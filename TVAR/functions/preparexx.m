function [ Y,X,Ystar ] = preparexx( data,L,tard,tarvar,VarBench )
Y=data;



%take lags
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];

%compute threshold variable
[Ystar0,Lx]=transform(Y(:,tarvar),VarBench.transform);
Ystar=lag0(Ystar0,tard);


Y=Y(max([L,tard,Lx+tard])+1:end,:);
X=X(max([L,tard,Lx+tard])+1:end,:);
Ystar=Ystar(max([L,tard,Lx+tard])+1:end,:);

end

