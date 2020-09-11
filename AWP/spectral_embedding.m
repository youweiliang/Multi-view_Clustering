function [indicator, eigvalue, L] = spectral_embedding(W, K)
% W = regularizedAdjacency(W, 2);
% d = (sum(W, 2));
% d = d ./ norm(d);
% d = d.^3;
% d = d ./ norm(d);
% D = sparse(1:size(W, 1), 1:size(W, 2), d);
% W = D*W*D;

degs = sum(W, 2);
D = sparse(1:size(W, 1), 1:size(W, 2), degs);
L = D - W;

degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
% calculate normalized Laplacian
L = D * L * D;
L = (L + L')/2;

[indicator, eigvalue] = eigs(L, K, 'sa');  % 'smallestreal'

end