function [C] = sparse_from_idx(elements, idx, n, m)
% Build a sparse matrix from elements according to idx
% idx is a logical matrix
% it put the elements onto C according to the idx
% size(C) == (n, m), length(elements) == nnz(idx)

[I, J, ~] = find(idx);
diag_idx = I==J;
elements(diag_idx) = 0;
C = sparse(I, J, elements, n, m);

end