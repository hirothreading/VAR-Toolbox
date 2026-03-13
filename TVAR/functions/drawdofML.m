function [ outnum ] = drawdofML( vg,vmean,vscale,omegainv,T,vprior )
%Alpha(theta,theta*)
   [ pnew ] = postv( vprior,vmean,omegainv,T) ;
   [ pold ] = postv( vprior,vg,omegainv,T) ;
    accept=min([exp(pnew-pold);1]);
   
M=vg;
V=vscale;
   %q(thetaG,theta*)
qterm=multivariatenormal(vmean,M,V);
outnum=log(accept)+(qterm);