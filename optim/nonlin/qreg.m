function f = qreg(R)
% Get a quadratic regularization function
%
%   f = qreg(R);   
%       returns a function handle, which represents the quadratic 
%       regularization function as follows. 
%
%       Case 1: If R is a scalar, then 
%
%           f(x) = (R/2) * ||x||^2
%
%       Case 2: If R is a vector of length d > 1, then
%
%           f(x) = \sum_i (R_i/2) * x_i^2
%
%       Case 3: If R is a positive definite matrix of size d x d, then
%
%           f(x) = (1/2) * x' * R * x
%
%       In each case, x can be a vector of length d, or multiple times of 
%       d. Suppose the length of x is m times d. Then it will partition
%       x into m sections, and sum the values yielded for each section.
%
%       As an objective funtion, x can be used in either of the following
%       way:
%
%           v = f(x);
%           [v, g] = f(x);
%           [v, g, H] = f(x);
%
%       Here, v, g, and H are respectively the function value, gradient
%       vector, and Hessian matrix evaluated at x. This syntax makes it
%       convenient to be incorporated in an optimization system.
%       (See the help of combfun).
%       

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% verify input

if ~(isfloat(R) && ndims(R) == 2 && isreal(R))
    error('qreg:invalidarg', 'R should be a numeric matrix');
end

%% main

if isscalar(R) 
    f = @(x) qregfun_s(x, R);
    
elseif isvector(R)    
    if size(R, 2) > 1
        R = R.';
    end
    f = @(x) qregfun_v(x, R);
    
elseif size(R,1) == size(R,2) 
    f = @(x) qregfun_m(x, R);

else
    error('qreg:invalidarg', 'The size of R is invalid.');
end


%% objective functions

function [v, g, H] = qregfun_s(x, R)

v = (x'*x) * (R/2);
if nargout >= 2
    g = R * x;
end
if nargout >= 3
    d = size(x,1);
    H = diag(constmat(1,d,R));
end

function [v, g, H] = qregfun_v(x, R)

d = size(R,1);
dx = size(x,1);

if dx == d
    g = R .* x;
    v = (x' * g) / 2;
    if nargout >= 3
        H = diag(R);
    end
    
elseif rem(dx, d) == 0
    m = dx / d;
    X = reshape(x, d, m);
    
    G = bsxfun(@times, R, X);
    v = sum(dot(X, G, 1)) / 2;
    if nargout >= 2
        g = G(:);
    end
    if nargout >= 3
        H = diag(repmat(R, m, 1));
    end    

else
    error('qreg:rterror', 'Dimension mismatch in qreg function.');
    
end
    
    

function [v, g, H] = qregfun_m(x, R)

d = size(R,1);
dx = size(x,1);

if dx == d
    g = R * x;
    v = (x' * g) / 2;
    if nargout >= 3
        H = R;
    end
    
elseif rem(dx, d) == 0
    m = dx / d;
    X = reshape(x, d, m);
    
    G = R * X;
    v = sum(dot(X, G, 1)) / 2;
    if nargout >= 2
        g = G(:);
    end
    if nargout >= 3
        H = zeros(dx, dx);
        for i = 1 : m
            sp = (i-1) * d + 1;
            ep = i * d;
            H(sp:ep, sp:ep) = R;
        end
    end
    
else
    error('qreg:rterror', 'Dimension mismatch in qreg function.');
    
end
    

