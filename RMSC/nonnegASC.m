function [ X ] = nonnegASC( B, opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% min |X-B|_F^2 
% s.t Xe=e, X>=0 where e is the constant one vector.

[n,m]=size(B);
A=repmat(1:m,n,1);
B_sort=sort(B,2,'descend');
cum_B=cumsum(B_sort,2);
sigma=B_sort-(cum_B-1)./A;
tmp=sigma>0;
idx=sum(tmp,2);
tmp=B_sort-sigma;
sigma=diag(tmp(:,idx));
sigma=repmat(sigma,1,m);
X=max(B-sigma,0);

end

