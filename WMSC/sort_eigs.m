function [eigvector, eigvalue] = sort_eigs(eigvector, eigvalue)
[D,I] = sort(diag(eigvalue), 'descend');
eigvalue = diag(D);
eigvector = eigvector(:, I);