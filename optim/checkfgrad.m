function checkfgrad(f, X, dx, op)
% A simple function for verifying that the gradient calculation is correct
%
%   checkfgrad(f, X, dx);
%       computes the objective values and gradients at the points given
%       by columns of X, and compares the evaluated gradient with the
%       the approximated gradient derived by finite difference.
%
%       Here, f is a function that supports the following syntax:
%
%           [v, g] = f(x);
%
%       Here, v and g are respectively the objective value and gradient
%       of f at x.
%
%       dx is the difference used in computing approximate gradient.
%
%   checkfgrad(f, X, dx, 'hess');
%       checks the calculation of Hessian matrix additionally. In this
%       case, the function f should support the following synax:
%
%           [v, g, H] = f(x);
%
%       Here, H is the Hessian matrix of f at x.
%

% Created by Dahua Lin, on Aug 7, 2010
%

%% main

check_hess = nargin >= 4 && strcmpi(op, 'hess');

n = size(X, 2);

for i = 1 : n
    
    x = X(:, i);
        
    if ~check_hess    
        g_a = approx_grad(f, x, dx);
        [v, g_e] = f(x); %#ok<ASGLU>    
        fprintf('at point %d:  norm(grad_diff) is %g\n', i, norm(g_a - g_e));    
    else
        g_a = approx_grad(f, x, dx);
        H_a = approx_Hessian(f, x, dx);
        
        [v, g_e, H_e] = f(x); %#ok<ASGLU> 
        fprintf('at point %d:  norm(grad_diff) is %g,  norm(hess_diff) is %g\n', ...
            i, norm(g_a - g_e), norm(H_a - H_e, 'fro'));
    end
end


%% the function to calculate approximate gradient

function g = approx_grad(f, x, dx)

d = length(x);
g = zeros(d, 1);

for k = 1 : d
    
    xp = x;
    xp(k) = xp(k) + dx;
    xn = x;
    xn(k) = xn(k) - dx;
    
    vp = f(xp);
    vn = f(xn);
    
    g(k) = (vp - vn) / (2 * dx);    
end


function H = approx_Hessian(f, x, dx)

d = length(x);
H = zeros(d, d);

for i = 1 : d
    for j = 1 : d
        if (i == j)
            xp2 = x; xp2(i) = xp2(i) + 2 * dx;
            xp1 = x; xp1(i) = xp1(i) + dx;
            xn1 = x; xn1(i) = xn1(i) - dx;
            xn2 = x; xn2(i) = xn2(i) - 2 * dx;
            
            H(i, i) = (-f(xp2) + 16 * f(xp1) - 30 * f(x) + 16 * f(xn1) - f(xn2)) / (12 * dx^2);
        else
            xpp = x; xpp(i) = xpp(i) + dx; xpp(j) = xpp(j) + dx;
            xpn = x; xpn(i) = xpn(i) + dx; xpn(j) = xpn(j) - dx;
            xnp = x; xnp(i) = xnp(i) - dx; xnp(j) = xnp(j) + dx;
            xnn = x; xnn(i) = xnn(i) - dx; xnn(j) = xnn(j) - dx;
            
            H(i, j) = (f(xpp) - f(xpn) - f(xnp) + f(xnn)) / (4 * dx^2);
        end
    end
end

