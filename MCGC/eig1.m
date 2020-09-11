function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)
warning('off', 'MATLAB:eigs:NotAllEigsConverged')
warning('off', 'MATLAB:eigs:IllConditionedA')
warning('off', 'MATLAB:eigs:SigmaNearExactEig')
if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end

if nargin < 3
    isMax = 1;
    isSym = 1;
end

if isMax
    error('Could not compute max eigenvalues')
end

if nargin < 4
    isSym = 1;
end

if isSym == 1
    A = max(A,A');
end

% if nnz(A)/numel(A) < 0.2
if nnz(A)/numel(A) < 0.1
	A = sparse(A);
elseif nnz(A)/numel(A) > 0.4
    A = full(A);
end
tol = 1e-3;
opts.tol = tol;
% try
    [v, d] = eigs(A, c + 2, eps, opts); % 'Tolerance',1e-3, 'MaxIter', 1000
% catch ME
%     warning(ME.message)
%     warning('Turn to using full eig function instead')
%     [v, d] = eig(full(A));
% end
not_conv_eigv = nnz(isnan(d));
cc = c + 2;
maxiter = 300;
while not_conv_eigv~=0
    tol = tol * 10;
    maxiter = maxiter * 2;
    cc = cc + c;
    [v, d] = eigs(A, cc, eps, 'Tolerance', tol, 'MaxIter', maxiter);
    found = nnz(diag(d) <= 0);
    if found >= c+1
        break
    end
    not_conv_eigv = nnz(isnan(d));
    if cc/c >= 5
        if size(A, 1) <= 5000
            [v, d] = eig(full(A));
        else
%             d(isnan(d)) = 0;
        end
        break
    end
%     fprintf('%d\n', not_conv_eigv)
end

% else
%     A = full(A);
%     [v, d] = eig(A);
% end

d = diag(d);

if isMax == 0
    [~, idx] = sort(d);
else
    [~, idx] = sort(d,'descend');
end
idx1 = idx(1:c);
eigval = d(idx1);
eigvec = v(:,idx1);
eigval_full = d(idx);
warning('on', 'MATLAB:eigs:NotAllEigsConverged')
warning('on', 'MATLAB:eigs:IllConditionedA')
warning('on', 'MATLAB:eigs:SigmaNearExactEig')