function [label] = WMSC_main(fea, numClust, knn0, b_hat, n_hat, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 6
    metric = 'squaredeuclidean';
end
WW = make_distance_matrix(fea, metric);
v = length(fea);
knn = knn0 + 1;

L = cell(1, v);
eigvalue = cell(1, v);
eigvector = cell(1, v);

opts.maxit = 2000;
opts.fail = 'keep';

for i = 1:v
    [A] = make_kNN_dist(WW{i}, knn);
    L{i} = normalizedLaplacian(A);  %% change from `L{i} = E - normalizedLaplacian(A);`
    [eigvector{i}, eigvalue{i}] = eigs(L{i}, numClust, 'la', opts);  %% change 'lm' to 'la'
    eigvalue{i}(isnan(eigvalue{i})) = 1;
    [eigvector{i}, eigvalue{i}] = sort_eigs(eigvector{i}, eigvalue{i});
end

label = WeightMSC(eigvector, eigvalue, L, v, numClust, b_hat, n_hat);
