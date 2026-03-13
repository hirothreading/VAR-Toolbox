function [Rmat,FindIndex]=kulian_jae_monthly_identification(SigmaMat,SignRes)

%BetaMat
ctiterion=1;
CholSigmaMat=chol(SigmaMat,'lower');
% iter=0;
FindIndex=1;
dyStar=sqrt(rows(SignRes));

while ctiterion > 0 %&& iter < 1000
    
    [~,Pmat]=var_sign_chol_restictions(CholSigmaMat(1:dyStar,1:dyStar),SignRes(:));
    Rmat=CholSigmaMat*blkdiag(Pmat,eye(rows(CholSigmaMat)-dyStar));
    
    SupElasticity1=Rmat(1,3)/Rmat(4,3);
    SupElasticity2=Rmat(1,2)/Rmat(4,2);
    
    DemElasticity=Rmat(1,1)/Rmat(4,1);
%     SuppIRF=varirf(BetaMat,Rmat(:,1),12);
    
    if SupElasticity1 < 0.0258 && SupElasticity2 < 0.0258 && -0.8<DemElasticity
        ctiterion=0;
    end
    
%     iter=iter+1;
end
% if iter==1000
%     FindIndex=0;
% end