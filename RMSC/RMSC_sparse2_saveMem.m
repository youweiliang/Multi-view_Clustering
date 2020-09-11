function [P] = RMSC_sparse2_saveMem(T, lambda, opts)
% We found that only input T can be sparse.
% Other matrices have high density during iterations!
% - Youwei Liang

%% OBJECTIVE
%min ||P||_* + lambda*\sum ||E_i||_1,1
% s.t P_i=P+E_i, Pe=e, P>=0 where e is the constant one vector.
%
% equal to ==>
%
% min ||Q||_* + lambda*\sum ||E_i||_1,1
% s.t P_i=P+E_i, Pe=e, P>=0, P=Q

%% RELATED PAPERS
% [1] Robust Multi-View Clustering via Low-rank and Sparse Decomposition.
%     Rongkai Xia, Yan Pan, Lei Du, and Jian Yin. In Proceedings of AAAI
%     Conference on Artificial Intelligence (AAAI), 2014

v = length(T);
[m, p] = size(T{1}); %num of samples
idx = logical(sparse(m, p));
for i=1:v
    idx = idx | sparse((T{i} ~= 0));
end

if m ~= p
    error('input matrix T must be a square matrix( transitions matrix ).\n');
end
Z = zeros(m, p);
E = cell(v, 1);
Y = cell(v, 1);
Q = zeros(m, p);
P = zeros(m, p);
for i=1:v
    E{i} = logical_extraction(randn(m, p), idx);
    Y{i} = zeros(m, p);
end

if isfield(opts, 'mu')
    mu = opts.mu;
else
    mu = 1e-3;
end
if isfield(opts, 'rho')
    rho = opts.rho;
else
    rho = 1.9;
end
if isfield(opts, 'max_iter')
    max_iter = opts.max_iter;
else
    max_iter = 100;
end

step = 0;
% P_old = randn(m, p);

while (1)
    % tic;
    step = step + 1;
    max_inf_norm = -1;
    for i = 1:v
        diff = T{i} - E{i} - P;
        inf_norm = norm(diff, 'inf');
        max_inf_norm = max(max_inf_norm, inf_norm);
    end
    %     funV=sum(svd(P))+lambda*norm(E(:),1);
%     relChg = norm(P-P_old, 'fro') / max(1, norm(P_old, 'fro'));
%     P_old = P;
    tmp = P - Q;
    max_inf_norm2 = norm(tmp(:), 'inf');
    if opts.DEBUG
        fprintf('iter %d: max_inf_norm = %s, relChg=%s, mu = %s, inf_norm2=%s, funV=%f\n', ...
            step, num2str(max_inf_norm), num2str(relChg), mu, num2str(max_inf_norm2), funV);
    end
    if step > 1 && max_inf_norm < opts.eps
        break;
    end

    if step > max_iter
        fprintf('reach max iterations %d \n', step);
        break;
    end

    % update P
    tmp = zeros(m, p);
    for i=1:v
        tmp = tmp + T{i}-E{i}-Y{i}/mu;
    end
    
    B = 1 / (v + 1) * (Q - Z / mu + tmp);
    P = proj_simplex(B);
    temp = sum(P, 2);
    if any(temp - 1.0 >= 1e-10)
        error('sum to 1 error');
    end

    % update Q
    M = P + Z / mu;
    C = 1 / mu;
    if svds(M, 1) < C
        Q = 0;
    else
        [U, Sigma, V] = svd(M, 'econ');
        Sigma = diag(Sigma);
        svp = length(find(Sigma > C));
        if svp >= 1
            Sigma = Sigma(1:svp) - C;
            Q = U(:, 1:svp) * diag(Sigma) * V(:, 1:svp)';
        else
            Q = 0;
        end
    end

    % update Ei
    for i = 1:v
        C = T{i} - P - Y{i} / mu;
        E{i} = max(C-lambda/mu, 0) + min(C+lambda/mu, 0);
%         E{i} = sparse_compare_max(C, lambda/mu) + sparse_compare_min(C, lambda/mu);
        Y{i} = Y{i} + mu * (P + E{i} - T{i});
    end
    Z = Z + mu * (P - Q);
    % update mu
    mu = min(rho*mu, 1e10);
end
if nnz(P)/numel(P) < 0.1
    P = sparse(P);
end
[pi, ~] = eigs(P', 1);
Dist = pi / sum(pi);
% pi = spdiags(Dist, 0, n, n);
pi05 = spdiags(Dist.^0.5, 0, m, m);
pi_05 = spdiags(Dist.^(-0.5), 0, m, m);

% P=(pi*P+P'*pi)/2;
P = pi05 * P * pi_05;
P = (P + P') / 2;
