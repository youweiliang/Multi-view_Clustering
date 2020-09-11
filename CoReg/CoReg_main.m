function [label] = CoReg_main(fea, numClust, knn0, lambda, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 5
    metric = 'squaredeuclidean';
end
numiter = 30;

WW = make_distance_matrix(fea, metric);
v = length(fea);
knn = knn0 + 1;

% Construct kernel
K = cell(v, 1);
for i = 1:v
    K{i} = make_kNN_dist(WW{i}, knn);
end

lambda = ones(v, 1) * lambda;
label = spectral_centroid_multiview_onkernel(K, v, numClust, lambda, numiter);
