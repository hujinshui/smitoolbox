function L = laplacemat(varargin)
% Compute the Laplacian matrix of a graph
%
%   L = laplacemat(g);
%   L = laplacemat(g, w);
%   L = laplacemat(g, w, a);
%
%       Computes the (regularized) Laplacian matrix for graph G. 
%
%       Suppose A is the affinity matrix of G, and d is the vector
%       of weighted degrees. Then the Laplacian matrix is defined
%       as follows:
%
%           L = diag(d + a) - A;
%
%       Input arguments:
%       - g:    the graph struct (yielded by make_gr), with g.dty = 'u'
%       - w:    the edge weights (can be either a scalar or a vector of
%               length m)
%       - a:    the regularization values to be added to the diagonal
%               entries (a scalar or a vector of length n)
%
%       In the output, L will be a sparse matrix of size n x n.
%
%   L = laplacemat(W);
%   L = laplacemat(W, a);
%
%       Computes the Laplacian matrix based on the edge weight matrix
%       W. Here, a is the regularization coefficient (scalar or
%       vector).
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 17, 2010
%       - Modified by Dahua Lin, on Nov 2, 2010
%           - support new graph structs.
%       - Modified by Dahua Lin, on Nov 13, 2010
%           - based on new graph class
%       - Modified by Dahua Lin, on Oct 27, 2011
%           - based on new graph struct.
%

%% verify input arguments

arg1 = varargin{1};

if isstruct(arg1)
    
    g = varargin{1};
    
    if nargin < 2
        w = 1;
    else
        w = varargin{2};
    end
        
    if nargin < 3
        a = [];
    else
        a = varargin{3};
    end
    
    if ~(is_gr(g) && g.dty == 'u')
        error('laplacemat:invalidarg', ...
            'g should be an undirected graph struct.');                
    end
    n = double(g.n);
    
    if ~(isfloat(w) && (isscalar(w) || (isvector(w) && numel(w) == g.m)))                
        error('laplacemat:invalidarg', ...
            'w should be a numeric vector of length g.m.');
    end
    
    is_wmat = false;

elseif isnumeric(arg1) || islogical(arg1)
    
    W = varargin{1};
    n = size(W, 1);
    if ~(ndims(W) == 2 && size(W, 2) == n)
        error('laplacemat:invalidarg', ...
            'W should be a square matrix.');
    end
    
    if nargin < 2
        a = [];
    else
        a = varargin{2};
    end    
    is_wmat = true;
    
else
    error('laplacemat:invalidarg', ...
        'The 1st argument to laplacemat is invalid.');
end

if ~isempty(a)
    if ~(isfloat(a) && (isscalar(a) || (isvector(a) && numel(a) == n)))
        error('laplacemat:invalidarg', ...
            'a should be a numeric vector of length n.');
    end
    if size(a, 2) > 1; a = a.'; end
    if ~isa(a, 'double'); a = double(a); end
end


%% main

if is_wmat    
    
    if issparse(W)
        if islogical(W)
            [i, j] = find(W);
            L = make_spL(n, i, j, 1, a);
        else
            [i, j, w] = find(W);
            L = make_spL(n, i, j, w, a);
        end
    else
        L = -W;
        dv = sum(W, 1);
        if ~isempty(a)
            dv = dv + a;
        end
        dind = 1 : n+1 : n*n;
        L(dind) = L(dind) + dv;
    end
    
else
    s = double(g.edges(1,:)).';
    t = double(g.edges(2,:)).';
    
    if isscalar(w)
        L = make_spL(n, s, t, w, a);
    else    
        if size(w, 2) > 1; w = w.'; end    
        L = make_spL(n, s, t, [w; w], a); 
    end
end


%% sub function

function L = make_spL(n, s, t, w, a)

if isscalar(w)
    dv = w * intcount(n, s).';
else
    dv = aggreg(w, n, s, 'sum');
end
if ~isempty(a)
    dv = dv + a;
end

m = numel(s);
if isscalar(w)
    w = w(ones(1,m), 1);
end
if isscalar(dv)
    dv = dv(ones(1,n), 1);
end

i = [s; (1:n).'];
j = [t; (1:n).'];
v = [-w; dv];

L = sparse(i, j, v, n, n);

