%% equivalent expressions may have different outcomes
% disp(eps);
format short
1+eps/2-1
1+(eps/2-1)

%%
format long    % show all the digits 
a = 1;  b = -(1e6+1e-6);  c = 1;
x1 = (-b + sqrt(b^2-4*a*c)) / (2*a)
x2 = (-b - sqrt(b^2-4*a*c)) / (2*a)
%%
% The first value is correct to all stored digits, 
% but the second has fewer % than six accurate digits: 
-log10(abs(1e-6-x2)/1e-6 )

%%
format long    % show all the digits
a = 1;  b = -(1e6+1e-6);  c = 1;
%%
% First, we find the ``good'' root using the quadratic forumla.  
x1 = (-b + sqrt(b^2-4*a*c)) / (2*a);
% Then we use the better formula for computing the other root.  
x2 = c/(a*x1)


