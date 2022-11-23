

%% ----------- Example 1 ------------------

%%
% We create two vectors for data about the population of China. The first
% has the years of census data, the other has the numbers of millions of
% people.
year = (1980:10:2010)'  
pop = [984.736; 1148.364; 1263.638; 1330.141];

%%
% It's convenient to measure time in years since 1980.
t = year - 1980;
y = pop;

%%
% Now we have four data points $(t_1,y_1),\ldots,(t_4,y_4)$, so $n=4$ and
% we seek an interpolating cubic polynomial. We construct the associated
% Vandermonde matrix: 
V = zeros(4,4);
for i = 1:4
    V(i,:) = [1 t(i) t(i)^2 t(i)^3];
end
V
 
%%
% To solve for the vector of polynomial coefficients, we use a backslash
c = V \ y

%%
% The algorithms used by the backslash operator are the main topic of this
% chapter. For now, observe that the coefficients of the cubic polynomial
% vary over several orders of magnitude, which is typical in this
% context. By our definitions, these coefficients are given in ascending
% order of power in $t$. MATLAB always expects the decreasing-degree order,
% so we convert ours to this convention here.
c = c(end:-1:1);       % reverse the ordering

%%
% We can use the resulting polynomial to estimate the population of China
% in 2005:
polyval(c,2005-1980)   % apply the 1980 time shift

%%
% The official figure is 1297.8, so our result is not bad. 

%%
% We can visualize the interpolation process. First, we 
% @glsbegin@plot@glsend@ the data as
% points. We'll shift the t variable back to actual years. 
plot(1980+t,y,'.')

%%
% We want to superimpose a plot of the polynomial. In order to add to a
% plot, we must use the @glsbegin@hold@glsend@ command:
hold on

%%
% To plot the interpolating polynomial, we create a vector with many points
% in the time interval using @glsbegin@linspace@glsend@.
tt = linspace(0,30,300)';   % 300 times from 1980 to 2010
yy = polyval(c,tt);         % evaluate the cubic
plot(1980+tt,yy)

%%
% Let's clear the figure and redo it, this time 
% continuing the curve outside of the original date range. We'll also 
% annotate the graph (using @glsbegin@title@glsend@, @glsbegin@xlabel@glsend@, 
% @glsbegin@ylabel@glsend@ and @glsbegin@legend@glsend@) to make its purpose 
% clear.
clf   % clear figure
plot(1980+t,y,'.')
hold on
tt = linspace(-10,50,300)';   
plot(1980+tt,polyval(c,tt))
title('Population of China')
xlabel('year'), ylabel('population (millions)')
legend('data','interpolant','location','northwest')

%%
% While the interpolation is plausible, the extrapolation to the future is
% highly questionable! As a rule, extrapolation more than a short distance
% beyond the original interval is not reliable.



%% ----------- Example 2 ------------------

%%
% Here is the system that ``broke" LU factorization for us.
A = [ 2 0 4 3; -2 0 2 -13 ; 1 15 2 -4.5 ; -4 5 -7 -10 ];
b = [ 4; 40; 29; 9 ];

%%
% When we use the built-in |lu| function with three outputs, we get the
% elements of the PLU factorization.
[L,U,P] = lu(A)

%%
% We can solve this as before by incorporating the permutation. 
x = backsub( U, forwardsub(L,P*b) )

%%
% However, if we use just two outputs with |lu|, we get
% $\mathbf{P}^T\mathbf{L}$ as the first result.
[PtL,U] = lu(A)

%%
% MATLAB has engineered the backslash so that systems with triangular _or
% permuted triangular_ structure are solved with the appropriate style of
% triangular substitution.
x = U \ (PtL\b)

%%
% The pivoted factorization and triangular substitutions are done silently
% and automatically when backslash is called on the original matrix. 
x = A\b


%% Example 3

%%
% MATLAB has a function cond to compute $\kappa_2(\bm{A})$. 
% The family of \index{matrix!Hilbert} _Hilbert matrices_ is famously 
% badly conditioned. Here is the
% $7\times 7$ case. 
A = hilb(7);
kappa = cond(A)

%%
% Next we engineer a linear system problem to which we know the exact answer.
x_exact = (1:7)';
b = A*x_exact;

%%
% Now we perturb the data randomly but with norm $10^{-12}$.
randn('state',333);          % reproducible results 
dA = randn(size(A));  dA = 1e-12*(dA/norm(dA));
% db = randn(size(b));  db = 1e-12*(db/norm(db));

%%
% We solve the perturbed problem using built-in pivoted LU and see how the
% solution was changed.
x = (A+dA) \ b; 
dx = x - x_exact;

%%
% Here is the relative error in the solution.
rel_error = norm(dx) / norm(x_exact)

%%
% And here are upper bounds predicted using the condition number of the
% original matrix. 
% b_bound = kappa * 1e-12/norm(b)
A_bound = kappa * 1e-12/norm(A)

%%
% Even if we don't make any manual perturbations to the data, machine
% epsilon does when we solve the linear system numerically.
x = A\b;
rel_error = norm(x - x_exact) / norm(x_exact)
rounding_bound = kappa*eps

%%
% Because $\kappa\approx 10^8$, it's possible to lose 8 digits of accuracy
% in the process of passing from $\bm{A}$ and $\bm{b}$ to
% $\bm{x}$. That's independent of the algorithm; it's inevitable once 
% the data are expressed in double precision. 

%%
% Now we choose an even more poorly conditioned matrix from this family.
A = hilb(14);
kappa = cond(A)

%%
% Before we compute the solution, note that $\kappa$ exceeds |1/eps|. In
% principle we might end up with an answer that is completely wrong.
rounding_bound = kappa*eps

%%
% MATLAB will notice the large condition number and warn us not to expect
% much from the result. 
x_exact = (1:14)';
b = A*x_exact;  x = A\b;

%%
% In fact the error does exceed 100%.
relative_error = norm(x_exact - x) / norm(x_exact)

