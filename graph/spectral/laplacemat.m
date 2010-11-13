function L = laplacemat(g, a)
% Compute the Laplacian matrix of a graph
%
%   L = laplacemat(G);
%   L = laplacemat(G, a);
%       Computes the (regularized) Laplacian matrix for graph G. 
%
%       Suppose A is the affinity matrix of G, and d is the vector
%       of weighted degrees. Then the Laplacian matrix is defined
%       as follows:
%
%           L = diag(d + a) - A;
%
%       Input arguments:
%       
%         - The graph G can be either an affinity matrix, or a object 
%           of class gr_edgelist (or its sub class).
%
%         - The regularization coefficient a can be in either of the 
%           following forms:
%               - a vector of length n
%               - a scalar
%               - omited or empty, when a is 0.
%
%       The returned matrix L will be a sparse matrix when G is 
%       a sparse matrix or a gr_edgelist object, or a dense matrix 
%       when G is a dense matrix.
%
%   Remarks
%   -------
%       - The caller should ensure that G is an undirected graph
%         (when input as a graph struct) or a symmetric matrix
%         (when input as an affinity matrix).
%
%       - The Laplacian matrix plays a central role in spectral graph
%         theory.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 17, 2010
%       - Modified by Dahua Lin, on Nov 2, 2010
%           - support new graph structs.
%       - Modified by Dahua Lin, on Nov 13, 2010
%           - based on new graph class
%

%% verify input arguments

if isa(g, 'gr_edgelist') && g.dtype == 'u'
    W = g.to_amat();
elseif isnumeric(g)
    W = g;   
    if ~(isfloat(W) && isreal(W) && ndims(W) == 2 && size(W,1) == size(W,2))
        error('laplacemat:invalidarg', 'W should be a real symmetric matrix.');
    end 
else
    error('laplacemat:invalidarg', 'The first argument is invalid.');
end

n = size(W, 1);

if nargin < 2 || isempty(a)
    a = 0;
else
    if ~(isfloat(a) && (isscalar(a) || (isvector(a) && numel(a) == n)))
        error('laplacemat:invalidarg', ...
            'a should be a real scalar or a real vector of length n.');
    end
    
    if size(a, 1) > 1   % turn a into a row vector
        a = a.';
    end
end

%% main

dv = sum(W, 1);

if issparse(dv)
    dv = full(dv);
end

if ~isequal(a, 0)
    dv = dv + a;
end

if issparse(W)
    L = spdiag(dv) - W;
else
    di = 1 + (0:n-1) * (n+1);
    L = -W;
    L(di) = L(di) + dv;
end


