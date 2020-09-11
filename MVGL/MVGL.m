function label = MVGL(S, Sv, S0, n, num_views, numClust, islocal_2)
% suppress warning for Updara_w
warning('off','MATLAB:nearlySingularMatrix')

J_old = 1; J_new = 10; EPS = 1e-3;
iter = 0;
while abs((J_new - J_old)/J_old) > EPS
    iter = iter +1;
    [~, A, ~] = Updata_A(S0, numClust,islocal_2);
    alpha = Updata_w(A,S);
    for i = 1:n
        for v = 1:num_views        
            Sv(:,i,v) = alpha(v,i).*S(:,i,v);
        end
    end
    S0 = sum(Sv,3);
    if iter == 30
        break
    end
    J_old = J_new;
    J_new =  sum(sum((A - S0).^2));
%     O(iter) = J_new;
end
AA = (A+A')/2;
% G = graph(AA);
% bins = conncomp(G);
% label = bins';
if nnz(AA)/numel(AA) < 0.1 && size(AA, 1) >= 3000
    AA = sparse(AA);
end
[label] = SpectralClustering(AA, numClust, 3);
warning('on','MATLAB:nearlySingularMatrix')