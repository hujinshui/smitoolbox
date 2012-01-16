function objf = comb_lossfun(X, Y, K, w, lossfun)
%COMB_LOSSFUN Combined loss function 
%
%   Generally, the total loss (risk) on a given set of (weighted) samples
%   is given by 
%
%       sum_{i=1}^n w_i loss(theta' * x_i, y_i).
%
%   Here, theta is a coefficient vector, x_i is the i-th sample (which
%   is associated with a weight w_i). y_i is the expected response. 
%   A loss value can be calculated using the function loss, which
%   takes two arguments: the prediction theta' * x_i, and the expected
%   response y_i.
%
%   Here, vector/matrix dimensions:
%   - x_i:          a sample vector: d x 1
%                   When there are n samples, they are grouped in a matrix
%                   of size d x n.
%
%   - theta:        the parameter vector/matrix: d x K.
%
%   - y_i:          a response vector: q x 1
%                   When there are n samples, all y_i vectors are grouped
%                   in a matrix of size q x n.
%
%   Note that in general, d and q need not be equal.
%
%
%   objf = COMB_LOSSFUN(X, Y, K, w, lossfun);
%
%       Constructs an objective function handle that evaluates the
%       total loss at all given samples for each parameter theta.
%
%       Suppose there are n training pairs of samples and responses.
%
%       Input arguments:
%       - X:        the sample (feature) matrix, size: d x n.
%
%       - Y:        the response matrix, size: q x n.
%
%       - K:        the number of columns in the parameter matrix.
%
%       - w:        the sample weights, size: 1 x n.
%                   w can also be [], indicating all sample weights are 1.
%
%       - lossfun:  the function handle of the loss function. It will be
%                   invoked in the following way:
%
%                   [v, g] = lossfun(theta' * X, Y).
%
%                   Here, v should be a 1 x n row vector containing the 
%                   objective values, and g should be a K x n matrix
%                   containing the gradient w.r.t. the first argument.
%
%                   
%       Output arguments:
%       - objf:     an objective function handle which can be used as
%                   follows:
%
%                       v = objf(theta).
%                       [v, g] = objf(theta).
%
%                   Here, v is the objective value, and g is the gradient
%                   vector of length d x K.
%
%       Note that objf can serve as an input to optimization functions
%       such as fminunc. The solution will be a vector of length d x K,
%       which is a concatenation of all columns of theta.
%       

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%


%% main

if size(w, 2) > 1
    w = w.';
end

objf = @(theta) rlmin_obj(theta, X, Y, K, w, lossfun);

%% The objective function

function [v, g] = rlmin_obj(theta, X, Y, K, w, lossfun)

d = size(X, 1);
if K > 1
    theta = reshape(theta, d, K);
end

Z = theta' * X;

out_g = nargout >= 2;

if out_g
    [L, Gz] = lossfun(Z, Y);  % Gz: K x n
else
    L = lossfun(Z, Y);
end

if isempty(w)
    v = sum(L);
    if out_g
        g = X * Gz';
    end
else
    v = L * w;
    if out_g
        if K == 1
            g = X * (Gz' .* w);
        else
            g = X * bsxfun(@times, Gz', w);
        end
    end    
end

if out_g && K > 1
    g = g(:);
end



