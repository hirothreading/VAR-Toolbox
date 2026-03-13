function [ vold,naccept ] = drawdof( vin,vscale,omegainv,T,vprior )
naccept=0;
vold=vin;
vnew=vin+randn(1,1).*sqrt(vscale);
if vnew<=0;
    accept=0;
else
   [ pnew ] = postv( vprior,vnew,omegainv,T) ;
   [ pold ] = postv( vprior,vin,omegainv,T) ;
   accept=exp(pnew-pold);
   
end
u=rand(1,1);
   if accept>u
       vold=vnew;
       naccept=1;
   end
end

