function [ lik,states,index1 ] = particlefilter(y,x,varcoef,iamat,Qmat,Fmat,npart,T,N,NS,L,B00,P00)
    EX=1;
	Cqmat=sqrt(Qmat);
	lik=0;
	states=zeros(T,NS);
    
    b0=repmat(B00,npart,1)+randn(npart,N)*chol(P00);
	
    for i=1:T  %loop through time
	%draw particles
	part=zeros(npart,NS);
    part(:,1:N)=randn(npart,N);	
    
	dens=zeros(npart,1);
	
	bnew=zeros(npart,NS);
	
    for j=1:npart
        
	xtemp=zeros(1,(N*L+EX));
    xtemp(1,1:(N*L))=x(i,1:(N*L));
    xtemp(1,cols(xtemp))=1;
    	
    %draw states
	
	htemp=(b0(j,:)*Fmat')+part(j,:)*Cqmat;
	bnew(j,:)=htemp;
	
	%compute conditional likelihood
	res=y(i,:)-xtemp*varcoef;
	sigma=iamat*diag(exp(htemp(1:N)))*iamat';
    isigma=invpd(sigma);
	dsigma=det(sigma);
	res2=res*isigma*res';
	dens(j)=(1/sqrt(dsigma))*exp(-0.5*res2);
    end
sdens=sum(dens);
prob=dens./sdens;
index1=discretesample(prob,length(prob));
b0=bnew(index1,:);
liki=sdens/npart;
lik=lik+log(liki);

states(i,:)=sum(b0.*repmat(prob,1,NS))';
    end

end

