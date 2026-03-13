function x = discretesample(p, n)
%==========================================================================
% Draw n samples (with replacement) from a discrete distribution.
%==========================================================================
% x = discretesample(p, n)
% -------------------------------------------------------------------------
% INPUT
%   - p : (1 x K) probability vector summing to 1
%   - n : number of samples to draw
% -------------------------------------------------------------------------
% OUTPUT
%   - x : (1 x n) vector of sampled indices in {1, ..., K}
% -------------------------------------------------------------------------
% Created by Dahua Lin, Oct 2008. Included in VAR Toolbox for TVAR support.
% =========================================================================

K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

% Construct bins from cumulative probabilities
edges = [0, cumsum(p)];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% Draw uniform random values and bin them
rv = rand(1, n);
c = histc(rv, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% Extract sample indices
xv = find(c);
if numel(xv) == n
    x = xv;
else
    xc = c(xv);
    d = zeros(1, n);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    x = cumsum(d);
end

x = x(randperm(n));
