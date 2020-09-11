function sparse_mat = build_from_idx(A, idx, n, width)
% size(A) = size(idx) = (n, width)

rowidx = ones(width, n) .* (1:n);
sparse_mat = sparse(rowidx, idx', A', n, n);