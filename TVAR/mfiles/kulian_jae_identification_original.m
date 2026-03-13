function [Rmat,FindIndex,normval]=kulian_jae_identification_original(SigmaMat,SignRes)

%BetaMat
ctiterion=1;
CholSigmaMat=chol(SigmaMat,'lower');
% iter=0;
FindIndex=1;

while ctiterion > 0 %&& iter < 1000
    
    [~,Pmat]=var_sign_chol_restictions(CholSigmaMat(1:4,1:4),SignRes(:));
    Rmat=CholSigmaMat*blkdiag(Pmat,eye(rows(CholSigmaMat)-4));
    
    SupElasticity1=Rmat(2,3)/Rmat(4,3);
    SupElasticity2=Rmat(2,2)/Rmat(4,2);
    
    DemElasticity=Rmat(2,1)/Rmat(4,1);
%     SuppIRF=varirf(BetaMat,Rmat(:,1),12);
    
    if SupElasticity1 < 0.6 && SupElasticity2 < 0.6 && -0.8<DemElasticity  
        ctiterion=0;
        normval=norm([SupElasticity1;SupElasticity2],2);
    end
    
%     iter=iter+1;
end
% if iter==1000
%     FindIndex=0;
% end