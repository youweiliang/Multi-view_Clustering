function [A, idx] = make_kNN_dist(D, knn)
% sigma = median(D);
n = size(D, 1);

[A, idx] = mink_new(D, knn, 2, 'sorting', false);
A = A(:);

rowidx = repmat((1:n)', knn, 1);
idx = idx(:);
non_diag = rowidx ~= idx;
rowidx = rowidx(non_diag);
idx = idx(non_diag);
A = A(non_diag);

sigma = mean(sqrt(A));
A = exp(-A/(2*sigma^2));

A = sparse(rowidx, idx, A, n, n);

% A = max(A, A');
A = (A + A')/2;
% A = power(A .* A', 1/2);

end