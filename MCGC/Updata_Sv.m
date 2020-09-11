function A = Updata_Sv(distX, c, k, islocal)
if islocal ~= 1
    error('islocal ~= 1')
end

NITER = 30;
num = size(distX,2);

[distX1, idx] = mink_new(distX, k+2, 2, 'sorting', true);

di = distX1(:,2:k+2);
rr = (k * di(:, k+1) - sum(di(:, 1:k), 2) + eps);
rr(rr==0) = eps;
tmp = (di(:, k+1) - di) ./ rr;
% for i = 1:num
%     di = distX1(i,2:k+2);
%     id = idx(i,2:k+2);
%     A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
% end
A = build_from_idx(tmp, idx(:,2:k+2), size(tmp, 1), size(tmp, 2));
lambda = 1;
A0 = (A+A')/2;

D0 = diag(sum(A0));
L0 = D0 - A0;
[F] = eig1(L0, c, 0);
% if sum(evs(1:c+1)) < 0.00000000001
%     error('The original graph has more than %d connected we component', c);
% end;
if islocal == 1
    idxa0 = idx(:,2:k+1);
else
    idxa0 = 1:num;
end
for iter = 1:NITER
    distf = L2_distance_1(F',F');
%     [distf1, ~] = sort(distf,2);
    if islocal == 1
        dfi = extract_from_idx(distf, idxa0);
    else
        dfi = distf(:,idxa0);
    end
    ad = -dfi/2/lambda;
    tmp = zeros(num, size(idxa0, 2));
    for i=1:num
        tmp(i, :) = EProjSimplex_new(ad(i, :));
    end
    if islocal == 1
        A = build_from_idx(tmp, idxa0, size(idxa0, 1), size(idxa0, 2));
    else
        warning('islocal ~= 1, compution may be wrong!\n')
        A = tmp;
    end
    A = (A+A')/2;
    D = diag(sum(A));
    L = D-A;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0);
%     evs(:,iter+1) = ev;

    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > 0.00000000001
        lambda = lambda/2;
    elseif fn2 < 0.00000000001
        lambda = lambda*2;  F = F_old;
    else
        break;
    end
    if iter == NITER
        warning('reach max iter')
    end
end