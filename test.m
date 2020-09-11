addpath ./utils
addpath ./MinMaxSelection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For text dataset (features are word frequence),
% use cosine metric

load('./data/MSRC_v1.mat')
numClust = length(unique(gt));
knn0 = 10;
metric = 'squaredeuclidean';

% load('./data/BBCSport.mat')
% numClust = length(unique(gt));
% knn0 = 10;
% metric = 'cosine';

cd ./AASC
[label] = AASC_main(fea, numClust, knn0, metric);
score = getFourMetrics(label, gt) %#ok<*NASGU,*NOPTS>

cd ..
cd ./AWP
[label] = AWP_main(fea, numClust, knn0, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./CoReg
lambda = 1e-1;
[label] = CoReg_main(fea, numClust, knn0, lambda, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./MCGC
beta = 1e-1;
[label] = MCGC_main(fea, numClust, knn0, beta, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./MVGL
[label] = MVGL_main(fea, numClust, knn0, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./RMSC
lambda = 1e-1;
[label] = RMSC_main(fea, numClust, knn0, lambda, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./WMSC
b_hat = 1e-1;
n_hat = 1e-1;
[label] = WMSC_main(fea, numClust, knn0, b_hat, n_hat, metric);
score = getFourMetrics(label, gt)

cd ..
cd ./SC
[label] = SC_multiview(fea, numClust, knn0, metric);
v = length(fea);
for i=1:v
    score = getFourMetrics(label{i}, gt)
end