function [B, D, rowidx, colidx]  = extract_from_idx(A, idx_matrix)
% Pick elements from A accroding to idx_matrix.
% Inputs:
%   A - the sourse matrix
%   idx_matrix - an index matrix, pick the elements from A if its colomn index in
%       the corespongding row of idx_matrix
% Outputs:
%   B - the extracted matrix, same shape as idx_matrix
%   D - the extracted matrix, same shape as A


[m, n] = size(idx_matrix);
rowidx = reshape(ones(n, m) .* (1:m), 1, n*m);
colidx = reshape(idx_matrix', 1, n*m);
idx = sub2ind(size(A), rowidx, colidx);
C = A(idx);
B = reshape(C, n, m);
B = B';

[i, j] = size(A);
D = sparse(rowidx, colidx, C, i, j);

% A = max(A, A');
% squareB = (D + D')/2;

end