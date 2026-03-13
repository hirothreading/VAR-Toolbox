function [irf0mat,Pmat]= var_sign_chol_restictions(cholsigmamat,signrest)

checka = 0;

dy     = size(cholsigmamat,1);
dimrest = size(signrest,1)/dy;

signrestmat = zeros(dy,dimrest);

for i = 1 : dimrest;
    signrestmat(:,i) = signrest((i-1)*dy+1:i*dy,:);
end
orderindices = zeros(1,dimrest);
signindices = zeros(dimrest,dimrest);
while checka < 1;
    pmat0 = randn(dy,dy);
    pmat = getqr(pmat0);
    terma = cholsigmamat*pmat;
    termaa = terma;
    TermA = 0;
    orderindices = zeros(1,dimrest);
    for i = 1 : dimrest
        checkb = 1;
        for j = 1 : dy;
            if isequal(sum(isnan(terma(:,j))),0)
                termb = terma(:,j) .* signrestmat(:,i) < 0;
                termc = -terma(:,j) .* signrestmat(:,i) < 0;
                termd = sum(termb);
                terme = sum(termc);
                if termd == 0 || terme == 0;
                    if termd == 0
                        checkb = 0;
                        orderindices(:,i) = j;
                        terma(:,j) = NaN*ones(dy,1);
                        signindices(:,i)=ones(dy,1);
                        break
                    else
                        checkb = 0;
                        terma(:,j) = NaN*ones(dy,1);
                        termaa(:,j) = -termaa(:,j);
                        orderindices(:,i) = j;
                        signindices(:,i)=-ones(dy,1);
                        break
                    end
                end;
            end
        end
        if isequal(checkb,1);
            break
        else
            TermA = TermA + 1;
        end
    end
    if isequal(TermA,dimrest);
        checka = 1;
    end
end

if ~isequal(dy,dimrest)
    TermB = 1:dy;
    TermC = zeros(1,dy-dimrest);
    TermD = 0;
    for i = 1 : dy
        if isempty(find(orderindices==TermB(:,i)))
            TermD = TermD + 1;
            TermC(:,TermD) = TermB(:,i);
        end
    end
    FinalOrder=[orderindices TermC];
else
    FinalOrder = orderindices;
end
irf0mat=termaa(:,FinalOrder);
Pmat = signindices.*pmat(:,FinalOrder);