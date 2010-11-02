function L = laplacemat(g, a)
% Compute the Laplacian matrix of a graph
%
%   L = laplacemat(g);
%   L = laplacemat(g, a);
%       computes the Laplacian matrix, which is defined as follows:
%
%           L = diag(d + a) - W.
%
%       Here, W the affinity matrix of the graph.
%
%       In the input, g can be a gr_edgelist (including gr_adjlist) 
%       struct, or an affinity matrix; a can be in either of the following 
%       form:
%       - a vector of length n
%       - a scalar
%       - omited or empty, when a is 0.
%
%       The returned matrix L will be a sparse matrix when g is 
%       a sparse matrix or a graph struct, or a dense matrix when
%       g is a dense matrix.
%
%   Remarks
%   -------
%       - The caller should ensure that g is an undirected graph
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
%

%% verify input arguments

if isstruct(g)
    W = edgelist_to_adjmat(g, 'u');    
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


