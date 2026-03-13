function out = iwpQ(v, ixpx)
%==========================================================================
% Draw from an Inverse-Wishart distribution.
%==========================================================================
% out = iwpQ(v, ixpx)
% -------------------------------------------------------------------------
% INPUT
%   - v    : degrees of freedom (scalar)
%   - ixpx : inverse of the scale matrix (k x k)
% -------------------------------------------------------------------------
% OUTPUT
%   - out  : (k x k) draw from IW(v, ixpx^{-1})
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit
% =========================================================================

k = size(ixpx, 1);
z = zeros(v, k);
for i = 1:v
    z(i,:) = (chol(ixpx)' * randn(k,1))';
end
out = inv(z' * z);
