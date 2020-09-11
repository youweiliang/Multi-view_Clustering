function [P] = RMSC(T, lambda, opts)

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

[m, p, n] = size(T); %num of samples
if m ~= p
    error('input matrix T must be a square matrix( transitions matrix ).\n');
end
Z = zeros(m, p);
E = randn(m, p, n);
Y = zeros(m, p, n);
Q = zeros(m, p);
P = zeros(m, p);

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
P_old = randn(m, p);

while (1)
    % tic;
    step = step + 1;
    max_inf_norm = -1;
    for i = 1:n
        Ti = T(:, :, i);
        Ei = E(:, :, i);
        diff = Ti - Ei - P;
        inf_norm = norm(diff, 'inf');
        max_inf_norm = max(max_inf_norm, inf_norm);
    end
    %     funV=sum(svd(P))+lambda*norm(E(:),1);
    relChg = norm(P-P_old, 'fro') / max(1, norm(P_old, 'fro'));
    P_old = P;
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
    B = 1 / (n + 1) * (Q - Z / mu + sum(T-E-Y/mu, 3));
    ASC_opts.tol = 1e-5;
    ASC_opts.DEBUG = 0;
    P = nonnegASC(B, ASC_opts);
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
    for i = 1:n
        C = T(:, :, i) - P - Y(:, :, i) / mu;
        E(:, :, i) = max(C-lambda/mu, 0) + min(C+lambda/mu, 0);
        Y(:, :, i) = Y(:, :, i) + mu * (P + E(:, :, i) - T(:, :, i));
    end
    Z = Z + mu * (P - Q);
    % update mu
    mu = min(rho*mu, 1e10);
end

[pi, ~] = eigs(P', 1);
Dist = pi / sum(pi);
pi = diag(Dist);
% P=(pi*P+P'*pi)/2;
P = (pi^0.5 * P * pi^-0.5 + pi^-0.5 * P' * pi^0.5) / 2;
