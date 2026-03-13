function [ temp] = safeexp( pfload )
sfactor=max(pfload);
temp=exp(pfload-sfactor);


end