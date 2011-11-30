function x = laplacesm(G, w, a, y)
% Performs Laplacian smooth based on a Gaussian MRF
%
%   The problem formulation is to minimize the following objective
%   function with respect to X:
%
%         (1/2) * sum_e w_e ||X(:,e_i) - X(:,e_j)||^2
%       + (1/2) * sum_i (1/2) * a_i ||X(:,i) - Y(:,i)||^2.
%
%   x = laplacesm(g, w, a, y);
%       solves the problem above using (regularized) Laplacian matrix.
%       
%       Inputs:
%       - G:    the input graph, in form of either an object of class
%               gr_edgelist, or an affinity matrix
%       - w:    the edge weight vector.
%       - a:    the regularization coefficients, in form of either
%               a scalar or a vector of length n.
%       - y:    the observation, which can be either a column vector of 
%               length n, or multiple column vectors arranged into an 
%               n x K matrix. 
%
%       Output:
%       - x:    a matrix of the same size as y (n x K), x(:,k) corresponds 
%               to the solution based on y(:,k).
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 17, 2010
%       - Modified by Dahua Lin, on Nov 2, 2010
%           - based on new graph struct
%       - Modified by Dahua Lin, on Nov 13, 2010
%           - based on new graph class
%       - Modified by Dahua Lin, on Nov 19, 2011
%           - based on new graph struct
%


%% verify input

if ~(isfloat(y) && isreal(y) && ndims(y) == 2 && size(y,1) == G.n)
    error('laplacesm:invalidarg', 'y should be a real matrix with n rows.');
end

%% main

L = laplacemat(G, w, a);  % this will verify the validity of g and a

if size(a, 2) > 1; a = a.'; end     % turns a into a column vector
    
K = size(y, 2);
if K == 1 || isscalar(a)
    ay = a .* y;
else
    ay = bsxfun(@times, a, y);
end

x = L \ ay;

