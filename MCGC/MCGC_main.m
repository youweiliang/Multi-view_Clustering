function [label] = MCGC_main(fea, numClust, knn0, beta, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 5
    metric = 'squaredeuclidean';
end
numiter = 5;
WW = make_distance_matrix(fea, metric);
v = length(WW);
n = length(WW{1});

% islocal_1 should be 1, because it is about local embedding - 2020.4.26
islocal_1 = 1;

A = zeros(n,n,v);
for i = 1:v 
    A(:,:,i) = Updata_Sv(WW{i}, numClust, knn0, islocal_1);
end

[label] = obj_MVCC(A,v,numClust,beta,numiter);