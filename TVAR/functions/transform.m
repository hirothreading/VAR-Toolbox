function [out,L]=transform(x,id)
if id==1
    logx=x;
    L=4;
    logxL=lag0(logx,L);
    out=(logx-logxL);
elseif id==2

    L=4;
    out=movsum(x,L-1);
elseif id==0
    out=x;
    L=0;
    
   
end