function affinity = make_affinity_D17()

s = sprintf('C:/youwei/ConsistentGraphLearning-master/pdist/D%d_data/sparse_w20.mat', 17);
load(s, 'sparse_w')
s = sprintf('C:/youwei/ConsistentGraphLearning-master/pdist/D%d_data/sigma.mat', 17);
load(s, 'sqrt_sigma', 'sigma')
% s = sprintf('D:/__MV-ProgressiveEC/Data/D17_data.mat');
% load(s, 'gt');
% truth = gt;
% numClust = 1000;
numView = 4;
% n_eig  = 1000;

affinity = cell(1, numView);

for v = 1:numView
    sig = sqrt_sigma(v);
    W = sparse_w{v};
    [i, j, s] = find(W);
    s = exp(-s /(2*sig^2));
    [m,n] = size(W);
    affinity{v} = sparse(i,j,s,m,n);
end
