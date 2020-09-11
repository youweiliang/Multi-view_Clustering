function [label] = WeightMSC(eigvector, eigvalue, L, v, c, b_hat, n_hat)
% if nargin == 3
%     b_hat = 0.1;
%     n_hat = 0.1;
% end
n = size(L{1}, 1);
Y = zeros(v);
T = cell(1, v);
largest_canonical_angle = zeros(v);

LV = cell(v,1);
for a=1:v
    T{a} = zeros(v);
    Vsigma = eigvector{a} * eigvalue{a};
    for i=1:v
        LV{i} = L{i} * eigvector{a};
        Y(i, a) = sum(sum(LV{i} .* Vsigma));
        largest_canonical_angle(i, a) = max(canonical_angle(eigvector{i}, eigvector{a}));
    end
    for i=1:v
        for j=1:v
            T{a}(i, j) = sum(sum(LV{i} .* LV{j}));
        end
    end
end

R = pi - largest_canonical_angle;
P = diag(sum(R, 2));
Q = P - R;
H1 = sum(cat(3,T{:}), 3);

beta = b_hat / sqrt(v) * norm(H1 + Q, 'fro');
eta = n_hat / norm(Q, 'fro') * norm(H1 + eye(v), 'fro') ;

H = H1 + beta*eye(size(H1)) + eta* Q;
H = (H+H')/2;
f = -sum(Y, 2)';
Aeq = ones(1, v);
beq = 1;
opts1 =  optimset('display','off');
x = quadprog(H,f,[],[],Aeq,beq,zeros(v,1),[], [], opts1);
% mu = reshape(x, 1,1,v);
% L_sum = cat(3, L{:}) .* mu;
% L_com = sum(L_sum, 3);
L_com = sparse(n, n);
for i=1:v
    L_com = L_com + x(i) * L{i};
end

[U, ~] = eigs(L_com, c, 'la');  %% add 'la' -- largest real eigenvalues
U = U ./ sqrt(sum(U.^2, 2));
label = kmeans(U, c, 'Replicates',10, 'MaxIter',1000);

end