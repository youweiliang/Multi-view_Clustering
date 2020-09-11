function [D, rowidx, colidx]  = extract_from_idx_ONE(idx_matrix, i, j)
% Pick elements from A accroding to idx_matrix with
% A being all one matrix
% Original function is: extract_from_idx(A, idx_matrix)
% Inputs:
%   A - the sourse matrix
%   idx_matrix - an index matrix, pick the elements from A if its colomn index in
%       the corespongding row of idx_matrix
% Outputs:
%   % B - the extracted matrix, same shape as idx_matrix
%   D - the extracted matrix, same shape as A


[m, n] = size(idx_matrix);
rowidx = reshape(ones(n, m) .* (1:m), 1, n*m);
colidx = reshape(idx_matrix', 1, n*m);
non_diag = rowidx ~= colidx;
rowidx = rowidx(non_diag);
colidx = colidx(non_diag);
D = sparse(rowidx, colidx, 1, i, j);

% A = max(A, A');
% squareB = (D + D')/2;

end