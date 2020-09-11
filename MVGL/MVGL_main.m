function [label] = MVGL_main(fea, numClust, knn0, metric)
% fea - a cell of feature, size(fea) = [numSample, numFeature]
% numClust - number of desired clusters
% knn0 - number of k-nearest neighbors
% metric - the metric for computing distance, can be (squaredeuclidean,cosine,original)
if nargin < 4
    metric = 'squaredeuclidean';
end

WW = make_distance_matrix(fea, metric);
v = length(fea);
n = length(WW{1});

% THIS SETTING FOLLOWS THE AUTHORS' CODE figure_03.m
islocal_1 = 1;
islocal_2 = 1;

S = zeros(n,n,v);
for i = 1:v 
    S(:,:,i) = Updata_Sv(WW{i}, numClust, knn0, islocal_1);
end
Sv = S./v;
S0 = sum(Sv,3);

label = MVGL(S, Sv, S0, n, v, numClust, islocal_2);
