function [g, H] = approx_deriv(f, x, h)
% Compute approximated derivatives by finite difference
%
%   g = approx_deriv(f, x);
%   g = approx_deriv(f, x, h);
%       computes the (approximated) gradient of a function f at x.
%       Here, h is the difference at x-value used in the approximation.
%       If omitted, it uses 1e-5 for h.
%
%   [g, H] = approx_deriv(f, x);
%   [g, H] = approx_deriv(f, h);
%       additionally computes the (approximated) Hessian matrix at x.
%
%   Remarks
%   -------
%       - First order central difference approximation is used.
%
%       - This routine is designed for the purpose of verifying whether
%         the computation of the derivatives of a particular function
%         is correctly implemented, but not for actual computation.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% verify input arguments

if ~isa(f, 'function_handle')
    error('approx_deriv:invalidarg', 'f should be a function handle.');
end

if ~(isfloat(x) && ndims(x) == 2 && size(x, 2) == 1)
    error('approx_deriv:invalidarg', 'x should be a numeric column vector.');
end
d = size(x, 1);

if nargin < 3
    h = 1e-5;
else
    if ~(isfloat(h) && isscalar(h) && isreal(h) && h > 0)
        error('approx_deriv:invalidarg', 'h should be a positive real number.');
    end
end

%% main

% compute gradient

g = zeros(d, 1);
for i = 1 : d
    x_p = av(x, i, h/2);
    x_n = av(x, i, -h/2);
    
    g(i) = (f(x_p) - f(x_n)) / h;
end

% compute Hessian

if nargout < 2
    return;
end

H = zeros(d, d);
for i = 1 : d
    for j = 1 : i
        x_pp = av(x, i, h/2, j, h/2);
        x_pn = av(x, i, h/2, j, -h/2);
        x_np = av(x, i, -h/2, j, h/2);
        x_nn = av(x, i, -h/2, j, -h/2);
        
        H(i, j) = (f(x_pp) + f(x_nn) - f(x_pn) - f(x_np)) / (h^2);
        if j < i
            H(j, i) = H(i, j);
        end
    end
end

%% Auxiliary routine

function y = av(x, i, ei, j, ej)

y = x;
y(i) = y(i) + ei;
if nargin >= 5
    y(j) = y(j) + ej;
end


