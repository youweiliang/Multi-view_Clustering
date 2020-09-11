function [C, L, U, D] = SpectralClustering(W, k, Type)
%SPECTRALCLUSTERING Executes spectral clustering algorithm
%   Executes the spectral clustering algorithm defined by
%   Type on the adjacency matrix W and returns the k cluster
%   indicator vectors as columns in C.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%   'Type' - Defines the type of spectral clustering algorithm
%            that should be used. Choices are:
%      1 - Unnormalized
%      2 - Normalized according to Shi and Malik (2000)
%      3 - Normalized according to Ng, Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering", 
%     Statistics and Computing 17 (4), 2007
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% calculate degree matrix
% d = (sum(W, 2));
% d = d ./ norm(d);
% d = d.^2;
% d = d ./ norm(d);
% D = sparse(1:size(W, 1), 1:size(W, 2), d);
% W = D*W*D;

n = size(W, 1);
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
L = D - W;

% compute normalized Laplacian if needed
switch Type
    case 2
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate inverse of D
        D_ = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D_ * L;
    case 3
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate D^(-1/2)
        D_ = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D_ * L * D_;
end

% L may be almost symmetric due to error in numerical calculation, so we make it symmetric.
L = (L + L')/2;

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
% try
    [U, ~] = eigs(L, k, 'sa');% eigs(L, speye(n), k, 'smallestreal');
% catch
%     warning('Error occur in eigs, the clustering result is wrong!\n')
%     U = rand(size(L,1), k);
% end
% eigcut(U, data)

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U_n = U ./ sqrt(sum(U.^2, 2));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
if any(isnan(U_n))
    U_n(isnan(U_n)) = 0;
end

if n < 5000
    C = kmeans(U_n, k,'Replicates',30,'MaxIter',1000);
elseif n < 10000
    C = kmeans(U_n, k,'Replicates',20,'MaxIter',1000);
elseif n < 30000
    C = kmeans(U_n, k,'Replicates',5,'MaxIter',1000);
else
    C = kmeans(U_n, k,'Replicates',3,'MaxIter',1000);
end
% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
% I = sparse(1:size(D, 1), C, 1);
% s = sqrt(sum(D * I, 1));
% H = I ./ s;

end