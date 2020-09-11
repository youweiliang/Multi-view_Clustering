function [label] = SC_multiview(fea, numClust, knn0, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 4
    metric = 'squaredeuclidean';
end

WW = make_distance_matrix(fea, metric);
v = length(fea);
knn = knn0 + 1;

label = cell(v, 1);

for i=1:v
    affn = make_kNN_dist(WW{i}, knn);
    label{i} = SpectralClustering(affn, numClust, 3);
end
