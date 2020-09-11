function [C] = logical_extraction(A, id)
% Extract elements from A according to the id
% id is a logical matrix, A is the target matrix,
% it extracts the elements in A corresponding to the nonzero elements in id
% size(C) == size(A) == size(id)

extracted = A(id);
linear_idx = find(id);
[n, m] = size(A);
B = sparse(linear_idx, ones(length(linear_idx),1), extracted, n*m, 1);
C = reshape(B, n, m);

end

