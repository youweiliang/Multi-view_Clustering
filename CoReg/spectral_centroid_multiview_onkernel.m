function [label] = spectral_centroid_multiview_onkernel(kernels, num_views, numClust, lambda, numiter)
tol = 1e-6;
[N , ~] = size(kernels{1});
opts.disp = 0;

numEV = numClust;
numVects = numClust;

L = cell(num_views, 1);
U = cell(num_views, 1);

objval = zeros(num_views+1, numiter+1);

for i=1:num_views
% Laplacian for the first view of the data
%    fprintf('computing kernel for X(%d)\n',i);
%    K(:,:,i) = constructKernel(X{i},X{i},options(i));
    %K1 = X1*X1';
    D = sum(kernels{i},1);
    D(D==0) = eps;
    %L1 = D1 - K1; 
    D = D.^-0.5;
    L{i} = D .* kernels{i} .* D';  
    L{i} = (L{i} + L{i}') / 2;
    [U{i}, E] = eigs(L{i}, numEV, 'LA', opts);    
    objval(i,1) = sum(diag(E));
end
objval(num_views+1, 1) = sum(objval(1:num_views, 1));

i = 2;
% now iteratively solve for all U's
while(i<=numiter+1)
    
    L_ustar = zeros(N);
    for j=1:num_views
        L_ustar = L_ustar + lambda(j)*U{j}*U{j}';
    end
    L_ustar = (L_ustar+L_ustar')/2;
    [Ustar, ~] = eigs(L_ustar, numEV, 'LA', opts);    

    L_ustar = Ustar*Ustar';
    L_ustar = (L_ustar+L_ustar')/2;
    for j=1:num_views            
        [U{j}, E] = eigs(L{j} + lambda(j)*L_ustar, numEV,'LA',opts);    
        objval(j,i) = sum(diag(E));
    end
    
    objval(num_views+1, i) = sum(objval(1:num_views, i));
    if i > 3 && abs(objval(num_views+1, i) - objval(num_views+1, i-1)) < tol
        break
    end
    
    i = i+1;
end

if (1)  %use view 1 in actual clustering
    U1 = Ustar;
    U1(isnan(U1)) = 0;
    normvect = sqrt(sum(U1.^2, 2));    

    U1 = U1 ./ normvect;
    U1(isnan(U1)) = 0;

    label = kmeans(U1(:,1:numVects),numClust,'Replicates',10,'MaxIter',1000);
end 

end