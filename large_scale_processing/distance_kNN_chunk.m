function [A, idx, A0] = distance_kNN_chunk(fea, metric, knn, n_chunk)
n = size(fea, 1);
batch_size = ceil(n / n_chunk);
batch_A = cell(1, n_chunk);
batch_idx = cell(1, n_chunk);
for i=1:n_chunk
    first_idx = (i-1)*batch_size + 1;
    last_idx = min(i*batch_size, n);
    dist_tmp = pdist2(fea(first_idx:last_idx, :), fea, metric);
    [batch_A{i}, batch_idx{i}] = mink_new(dist_tmp, knn, 2, 'sorting', false);
end
A = cat(1, batch_A{:});
idx = cat(1, batch_idx{:});
A0 = A;

% adjacency_matrix = zeros(n,n);
% for i=1:n
%     adjacency_matrix(i, idx(i,:)) = A(i, :);
% end
% A = sparse(adjacency_matrix);

rowidx = ones(knn, n) .* (1:n);
A = sparse(rowidx, idx', A', n, n);

% A = max(A, A');
A = (A + A')/2;
% A = power(A .* A', 1/2);
end