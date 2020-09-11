function [label] = obj_MVCC(W,num_views,numClust,beta,numiter)

N = size(W,1);
opts.disp = 0;
U = zeros(N,numClust,num_views);
gamma = 1;
S(1:N,1:N) = 0;
for v =1:num_views
%     fprintf('computing embedding matrix for view (%d)\n',v);
    [U(:,:,v),~] = eigs(W(:,:,v),numClust,'LA',opts);
    S = S + beta*U(:,:,v)*U(:,:,v)';
end
S = (S+S')/2;
DA = diag(sum(S));
LA = DA - S;
[H] = eig1(LA, numClust, 0);
zr = 10e-11;
k = 2;
OBJ = zeros(numiter+1,1);
for v =1:num_views
    OBJ(1) = OBJ(1) + trace(U(:,:,v)'*LA*U(:,:,v)) + norm(S-beta*U(:,:,v)*U(:,:,v)','fro');
end
while(k<=numiter+1)
%     fprintf('Running iteration %d\n',k-1);   
    A0 = zeros(N);
    for v = 1:num_views            
        [U(:,:,v), ~] = eigs(W(:,:,v) + beta.*S, numClust,'LA',opts);
        A0 = A0 + beta*U(:,:,v)*U(:,:,v)';
    end
    for iter = 1:50
        dist = L2_distance_1(H',H');
        S = A0.*0;
        for j = 1:N
            ai = A0(j,:);
            di = dist(j,:);
            ad = ai - 0.5.*gamma*di; 
            S(j,:) = EProjSimplex_new(ad);
        end
        S = (S + S.')/2;
        D = diag(sum(S));
        L = D - S;
        F_old = H;
        [H, ~, ev] = eig1(L, numClust, 0);
        fn1 = sum(ev(1:numClust));
        fn2 = sum(ev(1:numClust+1));
        if fn1 > zr
            gamma = gamma.*2;
        elseif fn2 < zr
            gamma = gamma/2;  H = F_old;
        else
            break;
        end
    end
    for v =1:num_views
        OBJ(k) = OBJ(k) + trace(U(:,:,v)'*L*U(:,:,v)) + norm(S-beta*U(:,:,v)*U(:,:,v)','fro');
    end
    if k > 3 && abs(OBJ(k) - OBJ(k-1)) < 1e-5
        break
    end
    if k == numiter
        fprintf('reach max iter in obj_MVCC\n')
    end
    k = k+1;
end
% plot(OBJ)
% [~, y]=graphconncomp(sparse(S)); 
% label = y';
if nnz(S)/numel(S) < 0.1 && size(S, 1) >= 500
    S = sparse(S);
end
[label] = SpectralClustering(S, numClust, 3);
