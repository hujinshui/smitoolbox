function f = mlogireg(X, y, w)
% Construct the objective for multinomial logistic regression
%
%   The objective function for an m-class logistic regression problem is
%
%       f(theta) = - sum_i w_i log p_i
%
%   Here, if the i-th sample is in the k-th class, then
%
%       p_i = exp(theta_k' * x_i) / sum_{l=1}^m exp(theta_l' * x_i)
%
%
%   f = mlogireg(X, y, w);
%       constructs the objective for multinomial logistic regression.
%
%       This function returns a function handle f which represent the
%       objective as formalized above.
%
%       Input arguments:
%       - X:        the design matrix. Suppose there are n input samples,
%                   each with q components, then X should be an matrix
%                   of size n x d.
%
%       - y:        y can be in either of the following forms.
%                   - a label vector of size n x 1. Here, y(i) can take
%                     value in {1, ..., m}. The function would set m
%                     to be max(y).
%                   - an assignment matrix of size n x m. Here, y(i,:)
%                     corresponds to X(i,:), and should satisfy the 
%                     properties: all(y(i,:) > 0) and sum(y(i,:)) = 1.
%
%       - w:        The weights of the samples. If all samples have the
%                   same weight, then w can be empty or omitted. 
%                   Otherwise, w should be a vector of length n. 
%
%       The output argument f is a function handle, which can be invoked 
%       as follows:
%
%           v = f(a);       
%           [v, g] = f(a);
%           [v, g, H] = f(a);
%
%       Here, f takes as input the parameter a, and returns the objective
%       value v, or optionally the gradient g, and Hessian matrix H.
%       This function handle can be used as an objective function in
%       numeric optimization.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 25, 2001.
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('mlogireg:invalidarg', 'X should be a numeric matrix.');
end
n = size(X, 1);

if ~(isnumeric(y) && ndims(y) == 2 && size(y, 1) == n)
    error('mlogireg:invalidarg', 'y should be a numeric matrix with n rows.');
end
m = size(y, 2);
if m == 1  % label vector
    m = max(y);
else
    if ~isfloat(y); y = double(y); end
end
   
if nargin < 3 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('mlogireg:invalidarg', 'w should be a vector of length n.');
    end
end

%% main

if size(y, 2) == 1   % use label vector
    f = genlr(@mlr_h, m, X, y, w);
    
else                 % use assignment matrix
    f = genlr(@mlr_ha, m, X, y, w);
    
end

%% h-functions

function [v, D1, D2] = mlr_h(U, y)
% use y as a label vector

n = size(U, 1);
P = nrmexp(U, 2);

si = (1:n).' + n * (y - 1); 
p = P(si);

v = - log(p);

if nargout >= 2
    D1 = P;
    D1(si) = D1(si) - 1;
end

if nargout >= 3
    D2 = make_H(P);
end


function [v, D1, D2] = mlr_ha(U, Y)
% use Y as an assignment matrix

P = nrmexp(U, 2);

v = - dot(log(P), Y, 2);

if nargout >= 2
    D1 = P - Y;
end

if nargout >= 3
    D2 = make_H(P);
end


%% Hessian matrix computation

function H = make_H(P)

[n, m] = size(P);
H = zeros(n, m, m);

for i = 1 : m
    pi = P(:, i);
    si = (i-1) * m + 1;
    ei = i * m;
    
    for j = 1 : i
        pj = P(:, j);
        sj = (j-1) * m + 1;
        ej = j * m;
                
        if i == j
            H(:,i,j) = (1 - pi) .* pi;
        else
            h = - pi .* pj;
            H(:,i,j) = h;
            H(:,j,i) = h;
        end
    end
end

    