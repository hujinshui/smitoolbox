function x = laplacesm(g, a, y)
% Performs Laplacian smooth based on a Gaussian MRF
%
%   The problem formulation is to minimize the following objective
%   function with respect to X:
%
%         (1/2) * sum_e w_e ||X(:,e_i) - X(:,e_j)||^2
%       + (1/2) * sum_i (1/2) * a_i ||X(:,i) - Y(:,i)||^2.
%
%   x = laplacesm(g, a, y);
%       solves the problem above using (regularized) Laplacian matrix.
%       In input, g is a weighted graph (in form of either graph struct
%       or affinity matrix), and a is a vector of regularization 
%       coefficients of length n. 
%
%       Y should can be a column vector of length n, or multiple
%       column vectors arranged into an n x K matrix. In the output,
%       X is a matrix of the same size, and X(:,k) corresponds to
%       the solution based on Y(:,k).
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 17, 2010
%       - Modified by Dahua Lin, on Nov 2, 2010
%           - based on new graph struct
%


%% main

L = laplacemat(g, a);  % this will verify the validity of g and a
if ~(isfloat(y) && ndims(y) == 2 && size(y,1) == g.n)
    error('laplacesm:invalidarg', 'The size of y is invalid.');
end

if size(a, 2) > 1; a = a.'; end     % turns a into a column vector
    
K = size(y, 2);
if K == 1 || isscalar(a)
    ay = a .* y;
else
    ay = bsxfun(@times, a, y);
end

x = L \ ay;

