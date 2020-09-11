function [y, S, evs] = Updata_A(A0, c, islocal)
if islocal ~= 1
    error('only deal with islocal == 1')
end
NITER = 60;
zr = 10e-11;
lambda = 1;

if nargin < 3
    islocal = 1;
end


A0 = A0-diag(diag(A0));
num = size(A0,1);
A10 = (A0+A0')/2;
D10 = diag(sum(A10));
L0 = D10 - A10;

% automatically determine the cluster number
[F0, ~, evs] = eig1(L0, c, 0);
% a = abs(evs); a(a<zr)=eps; ad=diff(a);
% ad1 = ad./a(2:end); 
% ad1(ad1>0.85)=1; ad1 = ad1+eps*(1:num-1)'; ad1(1)=0; ad1 = ad1(1:floor(0.9*end));
% [~, cs] = sort(ad1,'descend');
% sprintf('Suggested cluster number is: %d, %d, %d, %d, %d', cs(1),cs(2),cs(3),cs(4),cs(5))
% if nargin == 1
%     c = cs(1);
% end
F = F0(:,1:c);
% if sum(evs(1:c+1)) < zr
%     error('The original graph has more than %d connected component', c);
% end
if sum(evs(1:c)) < zr
    [~, y] = graphconncomp(sparse(A10)); y = y';
    S = A0;
    return;
end


for iter = 1:NITER
    dist = L2_distance_1(F',F');
    S = zeros(num);
    for i=1:num
        a0 = A0(i,:);
        if islocal == 1
            idxa0 = find(a0>0);
        else
            idxa0 = 1:num;
        end
        ai = a0(idxa0);
        di = dist(i,idxa0);
        ad = ai-0.5*lambda*di; 
        S(i,idxa0) = EProjSimplex_new(ad);
    end
    A = S;
    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0);
%     evs(:,iter+1) = ev;

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;  F = F_old;
    else
        break;
    end
    if iter == NITER
        warning('reach max iter in Updata_A')
    end
end

[clusternum, y]=graphconncomp(sparse(A)); y = y';
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', c)
end