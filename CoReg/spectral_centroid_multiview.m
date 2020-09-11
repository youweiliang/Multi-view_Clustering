function [UU] = spectral_centroid_multiview(X,num_views,numClust,lambda,truth,numiter, knn)
% INPUT:
% OUTPUT:

if (min(truth)==0)
    truth = truth + 1;
end

[N, M1] = size(X{1});
%[N M2] = size(X2);
K = zeros(N, N, num_views);
L = K;
U = zeros(N, numClust, num_views);

% for i=1:num_views
%     %options(i) = [];
%     options(i).KernelType = 'Gaussian';
%     options(i).t = sigma(i);
%     options(i).d = 4;
% end

opts.disp = 0;

numEV = numClust;
numVects = numClust;
for i=1:num_views
% Laplacian for the first view of the data
%     fprintf('computing kernel for X(%d)\n',i);
%     W = constructKernel(X{i},X{i},options(i), i);
    s = sprintf('kernel%d.mat', i);
    load(s);
    K(:,:,i) = full(kNN(W, knn));
    %K1 = X1*X1';
    D = diag(sum(K(:,:,i),1));
    %L1 = D1 - K1; 
    L(:,:,i) = sqrt(inv(D))*K(:,:,i)*sqrt(inv(D));  
    L(:,:,i)=(L(:,:,i)+L(:,:,i)')/2;
    [U(:,:,i), E] = eigs(L(:,:,i),numEV,'LA',opts);    
    objval(i,1) = sum(diag(E));
end

%%do clustering for first view
U1 = U(:,:,1);
normvect = sqrt(diag(U1*U1'));
normvect(normvect==0.0) = 1;
U1 = diag(normvect) \ U1;    

i = 2;
% now iteratively solve for all U's
while(i<=numiter+1)
%     fprintf('Running iteration %d\n',i-1);

    L_ustar(1:N,1:N) = 0;
    for j=1:num_views
        L_ustar = L_ustar + lambda*U(:,:,j)*U(:,:,j)';
    end
    L_ustar = (L_ustar+L_ustar')/2;
    [Ustar, Estar] = eigs(L_ustar, numEV,'LA',opts);    

    L_ustar = Ustar*Ustar';
    L_ustar = (L_ustar+L_ustar')/2;
    for j=1:num_views            
        [U(:,:,j), E] = eigs(L(:,:,j) + lambda*L_ustar, numEV,'LA',opts);    
        objval(j,i) = sum(diag(E));
    end

    objval(1,i) = sum(diag(E));

    if (1)  %use view 1 in actual clustering
        U1 = Ustar;
        normvect = sqrt(diag(U1*U1'));    
        normvect(normvect==0.0) = 1;
        U1 = diag(normvect) \ U1;
    end
    i = i+1;
end
UU = U1(:,1:numVects);

%%%CCA on U1 and U2
%i = i+1;
%[feats1 feats2 F_c P_c R_c nmi_c avgent_c] = multiviewccacluster(U1_norm, U2_norm, numClust, sigma1, sigma2, truth);
