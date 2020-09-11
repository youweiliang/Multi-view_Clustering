function [label] = AASC2020(WW, numClust, knn0)
% WW - full distance matrix
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
knn = knn0 + 1;
v = length(WW);
n = size(WW{1}, 1);

knn_idx = logical(sparse(n, n));
sigm = ones(v, 1);

for i=1:v
    W = WW{i};
%     sigm(i) = mean(sqrt(W), 'all');
    [~, idx] = mink_new(W, knn, 2);
    [tp] = extract_from_idx_ONE(idx, n, n);
    knn_idx = knn_idx | logical(tp);
end

aff = cell(1, v);
for i=1:v
    tmp_vec = WW{i}(knn_idx);
    sigm(i) = mean(sqrt(tmp_vec));
    tmp_vec = exp(-tmp_vec/(2*sigm(i)^2));
    tmp = sparse_from_idx(tmp_vec, knn_idx, n, n);
    aff{i} = (tmp + tmp') / 2; % need to be symmatric by requirement of AASC
end

label = AASC(aff,numClust);
