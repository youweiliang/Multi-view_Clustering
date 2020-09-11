function [label] = AWP_main(fea, numClust, knn0, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 4
    metric = 'squaredeuclidean';
end
WW = make_distance_matrix(fea, metric);
knn = knn0 + 1;
v = length(WW);

embedding = cell(1, v);
% clear WW W

opts.NITER = 300;
for t=1:v
    affn = make_kNN_dist(WW{t}, knn);
    embedding{t} = spectral_embedding(affn, numClust);
end
label = AWP(embedding, opts);