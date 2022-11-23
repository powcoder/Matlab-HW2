%%
% Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly
% dependent, we should find that the following matrix is somewhat
% ill-conditioned.
t = linspace(0,3,400)';
A = [ sin(t).^2, cos((1+1e-7)*t).^2, t.^0 ];
kappa = cond(A)

%%
% Now we set up an artificial linear least squares problem with a known
% exact solution that actually makes the residual zero.
x = [1;2;1];
b = A*x;  

%% 
% Using backslash to find the solution, we get a relative error that is
% about $\kappa$ times machine epsilon.
x_BS = A\b;
observed_err = norm(x_BS-x)/norm(x)
max_err = kappa*eps

%%
% If we formulate and solve via the normal equations, we get a much larger
% relative error. With $\kappa^2\approx 10^{14}$, we may not be left with
% more than about 2 accurate digits.
N = A'*A;
x_NE = N\(A'*b);
observed_err = norm(x_NE-x)/norm(x)
digits = -log10(observed_err)


%% Sparsity

% We examine the sparsity structure of matrix factorization on synthetic and real data.
% Data matrices are downloaded from SuiteSparse Matrix Collection
% http://faculty.cse.tamu.edu/davis/matrices.html

% prob = load('illc1850.mat');
% prob = load('well1033.mat');
% prob = load('well1850.mat');
% prob = load('Maragal_5.mat');
% A = prob.Problem.A;
% N = A' * A;

% 
A = sprand(5000, 1000, 0.002);
% 
N = A' * A;

nnz(A)
nnz(N)

t0 = cputime;
% L = chol(N,'lower');
[L, flag, P] = chol(N, 'lower');
t1 = cputime - t0;
% [L, D, p] = ldl(N);

fprintf('L nnz %.2f%%\n', 100*nnz(L)/(numel(L)/2 + m /2));

t0 = cputime;
[Q,R, P] = qr(A, 0);
t2 = cputime - t0;
fprintf('Q nnz %.2f%%\n', 100*nnz(Q)/numel(Q));
fprintf('R nnz %.2f%%\n', 100*nnz(R)/(numel(R)/2 + n / 2));
fprintf('nnz(Q+R) / nnz(L) %f\n', (nnz(Q) + nnz(R)) / nnz(L));
fprintf('t1 %.3e / t2 %.3e\n', t1, t2);

subplot(1,3,1); spy(L);
subplot(1,3,2); spy(Q);
subplot(1,3,3); spy(R);


