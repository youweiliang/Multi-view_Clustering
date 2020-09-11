function [X] = proj_simplex(B)

[n,m]=size(B);
B_sort=sort(B,2,'descend');
cum_B=cumsum(B_sort,2);
tmp = (1 - cum_B) ./ (1:m);
sig = B_sort + tmp;
sig = sig > 0;
[ir,ic] = find(sig);
maxcolind = accumarray(ir,ic,[n,1],@max);
% lambda = tmp(1:n, maxcolind);
idx = bsxfun(@eq, cumsum(ones(size(tmp)), 2), maxcolind);
lambda = sum(tmp.*idx, 2);
X = max(B + lambda,0);

end