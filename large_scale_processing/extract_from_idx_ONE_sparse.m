function [B, D, rowidx, colidx]  = extract_from_idx_ONE_sparse(sz_a, idx_matrix)
% Pick elements from A accroding to idx_matrix.
% Inputs:
%   A - the sourse matrix
%   idx_matrix - an index matrix, pick the elements from A if its colomn index in
%       the corespongding row of idx_matrix
% Outputs:
%   B - the extracted matrix, same shape as idx_matrix
%   D - the extracted matrix, same shape as A


[m, n] = size(idx_matrix);
rowidx = reshape(ones(n, m) .* [1:m], 1, n*m);
colidx = reshape(idx_matrix', 1, n*m);
idx = sub2ind([sz_a, sz_a], rowidx, colidx);
C = ones(length(idx), 1);
B = reshape(C, n, m);
B = B';

D = sparse(rowidx, colidx, C, sz_a, sz_a);

% A = max(A, A');
% squareB = (D + D')/2;

end