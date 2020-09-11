function U = baseline_spectral_onRW2(Prob,numClust,projev)
% spectral clustering using transition matrix K. Clustering is done by
% running k-means on the top-'numClust' eigenvectors of the normalized Laplacian
% INPUT:
% K: N x N similarity matrix. 
% numClust: desired number of clusters
% truth: N x 1 vector of ground truth clustering
% projev: number of top eigenvectors to return in V
% OUTPUT:
% V: top-'projev' eigenvectors of the Laplacian
% E: top-'projev' eigenvalues of the Laplacian
% F, P, R, nmi, avgent, AR: F-score, Precision, Recall, normalized mutual
% information, average entropy, Avg Rand Index
% C: obtained cluster labels 

% [p,~]=eigs(Prob',1);
% p=p/sum(p);
% P=(diag(p)*Prob+Prob'*diag(p))/2;
numEV = numClust*projev;
opts.disp = 0;
[V, E] = eigs(Prob,ceil(numEV),'lm',opts);  
U = V(:,1:ceil(numClust*1));

%[U E] = eig(L);   
%[E1 I] = sort(diag(E));  %sort in increasing order
%U = U(:,I(end-numEV+1:end));

norm_mat = repmat(sqrt(sum(U.*U,2)),1,size(U,2));
%%avoid divide by zero
for i=1:size(norm_mat,1)
    if (norm_mat(i,1)==0)
        norm_mat(i,:) = 1;
    end
end
U = U./norm_mat;


    