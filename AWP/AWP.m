function [Y, p, obj] = AWP(Fs, opts)
% solve:
% min_{R_i, Y} \sum_i || Y - F_i R_i ||_F
% s.t. R_i^T R_i = I, Y \in Ind
% 
% paras:
% Fs    1xv cell, contain nxk spectral embedding (F_i)
% opts  options for debug and parameter setting
% Y     nx1 indicator vector
% p vx1 vector
%
% Ref: Multiview Clustering via Adaptively Weighted Procrustes
% Written by: Lai Tian, tianlai.cs@gmail.com
% 2018-01-20

if (~exist('opts', 'var'))
    opts = [];
end

NITER = 20;
if (isfield(opts, 'NITER'))
    NITER = opts.NITER;
end
obj = zeros(NITER,1);
flag = 0;
obj_change=10;
v = length(Fs);
[n, k] = size(Fs{1});

% == initialization ==

%Y = sparse(1:n, randi(k, n, 1), ones(n, 1), n, k);
Y = zeros(n, k);
Rs = cell(v, 1);
for idx = 1:v; Rs{idx} = eye(k); end
p = ones(v, 1);

for iter = 1:NITER
    % == update Y ==
    sum_pFR = zeros(size(Y));
    for idx = 1:v
        sum_pFR = sum_pFR + Fs{idx}*Rs{idx}/p(idx);
    end
    [~, y_idx] = max(sum_pFR, [], 2);
    Y = full(sparse(1:n, y_idx, ones(n, 1), n, k));

    % == update R_i ==
    for idx = 1:v
       [tmp_u, ~, tmp_v] = svd(Y'*Fs{idx});
       Rs{idx} = tmp_v * tmp_u';
    end

    % == update p_i ==
    for idx = 1:v
       p(idx) = norm(Y - Fs{idx}*Rs{idx}, 'fro'); 
    end

    obj(iter) = norm(p);
    if iter > 1
        obj_change = min(abs(obj(iter)-obj(1:iter-1)))/abs(obj(1) - obj(iter));
    end
    if obj_change < 1e-6
        if flag == iter-1 % two consecutive obj_change < 1e-6
            break;
        end
        flag = iter;
    end
end

% == from indicator mat to vec ==

[~, Y] = max(Y, [], 2);

end