function out=dlog(x)
    logx=log(x);
  
    logxL=lag0(logx,1);
    out=(logx-logxL)*100;
