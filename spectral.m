%adjacency matrix A from Graph G
A=[ 0 0.8 0.6 0.1 0 0; 0.8 0 0.9 0 0 0;
0.6 0.9 0 0 0 0.2; 0.1 0 0 0 0.6 0.7;
0 0 0 0.6 0 0.8;
0 0 0.2 0.7 0.8 0];
imagesc(A),title('Adjacency Matrix')
colorbar

% eigenvector to corresponding eigenvalue
[v d]=eig(A);

% laplacian matrix for the graph
figure;

Lg= [ 1.5 -0.8 -0.6 -0.1 0 0;
-0.8 1.7 -0.9 0 0 0;
-0.6 -0.9 1.7 0 0 -0.2;
-0.1 0 0 1.4 -0.6 -0.7;
0 0 0 -0.6 1.4 -0.8;
0 0 -0.2 -0.7 -0.8 1.7];

% image of the laplacian
imagesc(Lg)
colorbar

subplot(1,2,1)
imagesc(A),title('Adjacency Matrix'),colorbar
subplot(1,2,2)
imagesc(Lg),title('Laplacian Matrix'),colorbar
%LU matrix factorization
[L,U,P] = lu(Lg);
%Calculate the LU factorization of A by calling lu with three outputs.
%Generate spy plots of the L and U factors.
subplot(1,2,1)
spy(L)
title('L factor')
subplot(1,2,2)
spy(U)
title('U factor')


[ev ed]=eig(Lg);
% 'lambda2' is the Algebric connectivity and eigenvector 'v2' for
% corresponding lambda2
v2 = ev(:,2);
disp(v2)
plot(v2, 'r*')
hold on
plot(v2(1:3,:),  'r')
hold on
plot(v2(4:6,:), 'b')

% plot the eigen vector corresponding to the 2nd smallest eigen value
figure,plot(v2,'r*'),title('Eigenvector corresponding to the 2nd Smallest Eigenvalue');

 %s1=sort(diag(ed),'descend');
 %plot(s1)
%#### normalized cut ####------------------------

D= [ 1.5 0 0 0 0 0;
0 1.7 0 0 0 0;
0 0 1.7 0 0 0;
0 0 0 1.4 0 0;
0 0 0 0 1.4 0;
0 0 0 0 0 1.7];
 imagesc(D)
 colorbar
 Lnor=(D^-0.5).*Lg.*(D^-0.5);
 imagesc(Lnor)
 colorbar
 subplot(1,2,1)
imagesc(D),title('Diagonal Matrix'),colorbar
subplot(1,2,2)
imagesc(Lnor),title('Normalized Laplacian Matrix'),colorbar
 
Lnorm= [ 1 -0.5 -0.4 -0.1 0 0;
-0.5 1 -0.5 0 0 0;
-0.4 -0.5 1 0 0 -0.1;
-0.1 0 0 1 -0.4 -0.5;
0 0 0 -0.4 1 -0.5;
0 0 -0.1 -0.5 -0.5 1];
imagesc(Lnorm),colorbar
title('clustering with normalized Laplacian')
[vn en]=eig(Lnorm)
%-----------------Condition number of matrix-------------
 
 cond(A) % check it is illied conditioned 
 
% spectral radius of the matrix 
%max(sqrt(eig(A'.*A)))
 
%% improving ill conditioned problem 
%adjacency matrix A from Graph G
A1=[  1.5944 0.8 0.6 0.1 0 0; 0.8  1.5944 0.9 0 0 0;
0.6 0.9   1.5944 0 0 0.2; 0.1 0 0   1.5944 0.6 0.7;
0 0 0 0.6   1.5944 0.8;
0 0 0.2 0.7 0.8   1.5944];
imagesc(A1),title('Adjacency Matrix')
colorbar

% eigenvector to corresponding eigenvalue
[v1 d1]=eig(A1);

% laplacian matrix for the graph
figure;

Lg= [ 1.5 -0.8 -0.6 -0.1 0 0;
-0.8 1.7 -0.9 0 0 0;
-0.6 -0.9 1.7 0 0 -0.2;
-0.1 0 0 1.4 -0.6 -0.7;
0 0 0 -0.6 1.4 -0.8;
0 0 -0.2 -0.7 -0.8 1.7];

% image of the laplacian
imagesc(Lg)
colorbar


%LU matrix factorization
[L,U,P] = lu(Lg);
%Calculate the LU factorization of A by calling lu with three outputs.
%Generate spy plots of the L and U factors.
subplot(1,2,1)
spy(L)
title('L factor')
subplot(1,2,2)
spy(U)
title('U factor')


[ev ed]=eig(Lg);
% 'lambda2' is the Algebric connectivity and eigenvector 'v2' for
% corresponding lambda2
v2 = ev(:,2);
disp(v2)
plot(v2, 'r*')
hold on
plot(v2(1:3,:),  'r')
hold on
plot(v2(4:6,:), 'b')

% plot the eigen vector corresponding to the 2nd smallest eigen value
figure,plot(v2,'r*'),title('Eigenvector corresponding to the 2nd Smallest Eigenvalue');

 %s1=sort(diag(ed),'descend');
 %plot(s1)
 
%Condition number of matrix
 
 cond(A) % check it is illied conditioned 
 
% spectral radius of the matrix 
max(sqrt(eig(A.'*A)))