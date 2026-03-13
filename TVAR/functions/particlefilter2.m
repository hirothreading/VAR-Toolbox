function [ lik,states,hstates ] = particlefilter2(y,varcoef,iamat,Qmat,Fmat,npart,T,N,L,EX,B00,P00,scle,b00,p00)
[F,MU]=comp(varcoef,N,L,EX);
NS=cols(F);
R=eye(N)*scle;
iR=invpd(R);
detR=det(R);
H=eye(N,NS);
Cqmat=chol(Qmat);
	lik=0;
	states=zeros(T,NS);
    hstates=zeros(T,N);
    
    h0=repmat(B00,npart,1)+randn(npart,N)*chol(P00);
    b0=repmat(b00,npart,1)+[randn(npart,N) zeros(npart,NS-N)]*chol(p00);
	
    for i=1:T  %loop through time
	%draw particles
	part=randn(npart,N);	
    hpart=randn(npart,N);
    
	dens=zeros(npart,1);
	
	bnew=zeros(npart,NS);
	hnew=zeros(npart,N);
    for j=1:npart
        
	
    %draw volatility states
	htemp=(h0(j,:)*Fmat')+hpart(j,:)*Cqmat;

	hnew(j,:)=htemp;
    %draw data states
    sigma=iamat*diag(exp(htemp(1:N)))*iamat';
    A0=cholx(sigma);
    reduced=part(j,:)*A0;

    errors=zeros(1,NS);
    errors(1,1:N)=reduced;
    bnew(j,:)=MU'+b0(j,:)*F'+errors;
	
	%compute conditional likelihood
	res=y(i,:)-(H*bnew(j,:)')';
        res2=res*iR*res';
	tempdens=(1/sqrt(detR))*exp(-0.5*res2);
	 if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens(j)=exp(-100);
disp(sprintf('warning= %s ', num2str([jj j])));

    else
        dens(j)=tempdens;
    end
    end
sdens=sum(dens);
prob=dens./sdens;
index1=discretesample(prob,length(prob));
b0=bnew(index1,:);
h0=hnew(index1,:);
liki=sdens/npart;
lik=lik+log(liki);

states(i,:)=sum(b0.*repmat(prob,1,NS))';
hstates(i,:)=sum(h0.*repmat(prob,1,N))';
    end

end