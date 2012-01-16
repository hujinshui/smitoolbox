function f = mlogiregf(X, y, w, rc)
%MLOGIREGF Multi-class logistic regression objective function
%
%   The objective function for an K-class logistic regression problem is
%
%       f(theta) = - sum_i w_i log p_i + (rc/2) * ||theta||^2
%
%   Here, if the i-th sample is in the k-th class, then
%
%       p_i = exp(z_{ki}) / sum_{l=1}^m exp(theta_l' * z_{li})
%
%   where z_{ki} = theta_k' * x_i + theta_{k0}
%
%
%   f = MLOGIREGF(X, y, w, rc);
%       constructs the objective for multi-class logistic regression.
%
%       This function returns a function handle f which represent the
%       objective as formalized above.
%
%       Input arguments:
%       - X:        The sample matrix of size d x n.
%
%       - y:        y can be in either of the following forms.
%                   - a label vector of size 1 x n. Here, y(i) can take
%                     value in {1, ..., m}. The function would set m
%                     to be max(y).
%                   - an assignment matrix of size K x n. Here, y(:,i)
%                     corresponds to X(:,i). Note that each column of y
%                     must sums to one.
%
%       - w:        The weights of the samples. If all samples have the
%                   same weight, then w can be empty or omitted. 
%                   Otherwise, w should be a vector of length n. 
%
%       - rc:       The regularization coefficient.
%
%       The output argument f is a function handle, which can be invoked 
%       as follows:
%
%           v = f(a);       
%           [v, g] = f(a);
%
%       Here, f takes as input the parameter a, and returns the objective
%       value v, or optionally the gradient g.
%       This function handle can be used as an objective function in
%       numeric optimization.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 25, 2001.
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('mlogiregf:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~(isfloat(y) && isreal(y) && ndims(y) == 2 && size(y, 2) == n)
    error('mlogiregf:invalidarg', 'y should be a real matrix with n columns.');
end
K = size(y, 1);
if K == 1  % label vector
    K = max(y);
end
   
if nargin < 3 || isempty(w)
    w = [];
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('mlogiregf:invalidarg', 'w should be a vector of length n.');
    end
end

if nargin < 4
    rc = 0;
else
    if ~(isfloat(rc) && isscalar(rc) && rc >= 0)
        error('mlogiregf:invalidarg', 'rc should be a non-negative scalar.');
    end
end


%% main

Xa = [X; ones(1, n)];
f_loss = comb_lossfun(Xa, y, K, w, @mlogistic_loss);

if rc == 0
    f = f_loss;
else
    c = [ones(d, K); zeros(1, K)];
    f_reg = tikregf(c(:));
    f = comb_objfun(1, f_loss, rc, f_reg);
end

    