function L = normalizedLaplacian(A)
% Fixed by Youwei. 2020.06.17

degs = sum(A, 2);
D = sparse(1:size(A, 1), 1:size(A, 2), degs);

degs(degs == 0) = eps;
% calculate D^(-1/2)
D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
% calculate normalized Laplacian
L = D * A * D;
L = (L + L')/2;