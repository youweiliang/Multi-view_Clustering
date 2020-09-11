function [label] = AASC_main(fea, numClust, knn0, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 4
    metric = 'squaredeuclidean';
end
WW = make_distance_matrix(fea, metric);

[label] = AASC2020(WW, numClust, knn0);

