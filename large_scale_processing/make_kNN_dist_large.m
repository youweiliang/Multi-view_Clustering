function make_kNN_dist_large(fea, save_path, metric)
addpath ../MinMaxSelection

mkdir(save_path)
v = length(fea);

n = length(gt);
knn_idx = sparse(n, n);
sparse_w = cell(1, v);
mink = cell(1, v);
for i=1:v
%     W  = pdist2(fea{i}, fea{i}, 'squaredeuclidean');
    [sparse_w{i}, idx, mink{i}] = distance_kNN_chunk(fea{i}, metric, K, 10);
    [~, tp] = extract_from_idx_ONE_sparse(n, idx);
    knn_idx = knn_idx + tp;
%     clear W
end

s = sprintf('%s/mink.mat', save_path);
save(s, 'mink')
s = sprintf('%s/sparse_w20.mat', save_path);
save(s, 'sparse_w')
knn_idx = logical(knn_idx + knn_idx');
s = sprintf('%s/knn_idx.mat', save_path);
save(s, 'knn_idx')
knn_distance = cell(1, v);
knn_dense = zeros(v, nnz(knn_idx));
for i=1:v
    W  = pdist2(fea{i}, fea{i}, 'squaredeuclidean');
    knn_distance{i} = logical_extraction(W, knn_idx);
    knn_dense(i, :) = W(knn_idx);
    clear W
end
s = sprintf('%s/knn_sparse.mat', save_path);
save(s, 'knn_distance')
s = sprintf('%s/knn_dense.mat', save_path);
save(s, 'knn_dense')
% ne = nnz(knn_idx);
% D = zeros(v, ne);
% for i=1:v
%     D(i,:) = knn_distance{i}(knn_idx);
% end
% s = sprintf('%s/D.mat', save_path);
% save(s, 'D')
