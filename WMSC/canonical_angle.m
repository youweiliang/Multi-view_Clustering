function zeta = canonical_angle(V1, V2)
A = V1' * V2;
[~,S,~] = svd(A);
zeta = acos(diag(S));