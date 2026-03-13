function out = safeexp(x)
%==========================================================================
% Numerically stable exponentiation: subtract the max before exp to avoid
% overflow.
%==========================================================================
% out = safeexp(x)
% -------------------------------------------------------------------------
% Adapted from Haroon Mumtaz's TVAR Toolkit
% =========================================================================

out = exp(x - max(x));
