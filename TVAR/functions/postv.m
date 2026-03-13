function [ pnew ] = postv( vprior,vdraw,lamda,N )

nu=(1/vprior)+0.5*sum((log(1./lamda)+lamda));
pnew=(N*vdraw/2)*log(0.5*vdraw)-N*gammaln(0.5*vdraw)+(-nu*vdraw);


end

