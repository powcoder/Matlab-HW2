
%% Bucky
% Create a 60 by 60 sparse matrix B and compute the power terms
[B,V] = bucky;
subplot(2,2,1), spy(B)
subplot(2,2,2), spy(B^2)
subplot(2,2,3), spy(B^2)
subplot(2,2,4), spy(B^3)


%% Cholesky factorization of sparse matrix

% Generate sparse matrix

n = 60;
A = zeros(n, n);

% Set diagonal as one
% A(1:1+size(A,1):end) = 1;
% A(:,1) = 1;
% A(1, :) = 1;
% subplot(1,3,1), spy(A)
% [L,U] = lufact(A);
% subplot(1,3,2), spy(L)
% subplot(1,3,3), spy(U)



A(1:1+size(A,1):end) = 1;
A(:,end) = 1;
A(end, :) = 1;
subplot(1,3,1), spy(A)
[L,U] = lufact(A);
subplot(1,3,2), spy(L)
subplot(1,3,3), spy(U)