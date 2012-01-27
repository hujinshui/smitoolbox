function edges = gr_prim(G, w, rv, msk, max_size, max_w)
%GR_PRIM Prim's Minimum Spanning Tree Algorithm
%
%   edges = GR_PRIM(G, w, rv);
%   edges = GR_PRIM(G, w, rv, msk);
%   edges = GR_PRIM(G, w, rv, msk, max_size);
%   edges = GR_PRIM(G, w, rv, msk, max_size, max_w);
%
%       Use Prim's algorithm to find a minimum spanning tree with a
%       specified root vertex over the given graph G.
%
%       The Prim's algorithm grows the tree from a single vertex rv,
%       gradually incorporating the edges that connect between a vertex
%       that has been in a tree and one that out of the tree.
%
%       Input arguments:
%       - G:        The graph, which should be a undirected graph with
%                   neighbor system.
%
%       - w:        The edge weights (a vector of length G.m).
%
%       - rv:       The root vertex
%
%       - msk:      The vertex mask. The vertex v is allowed, when
%                   msk(v) is true. If msk is not specified, then by
%                   default, all vertices are allowed.
%
%       - max_size: The maximum number of vertices in the tree.
%                   (If not specified, the function grows the tree to
%                    span the entire connected component that contains rv)
%
%       - max_w:    The maximum weight of edges that can be incorporated
%                   in the tree. (If not specified, all edges are allowed)
%
%       Note that msk, max_size, and max_w are optional parameters. If
%       don't want to enforce one of such constraints, you can input
%       the corresponding parameter as an empty array [].
%
%       Output arguments:
%       - edges:    The vertor comprised of the indices of the edges
%                   included in the tree (in nondecreasing order of 
%                   edge weights)
%

% Created by Dahua Lin, on Jan 26, 2012
%


%% verify input arguments

if ~(is_gr(G) && isequal(G.dty, 'u') && G.has_nbs)
    error('gr_prim:invalidarg', ...
        'G should be an undirected graph struct (with neighbor system).');
end
n = G.n;

if ~(isfloat(w) && isreal(w) && ~issparse(w) && numel(w) == G.m)
    error('gr_prim:invalidarg', ...
        'w should be a real vector of length G.m.');
end
w = [w w];

if ~(is_pos_scalar(rv) && rv == fix(rv) && rv <= n)
    error('gr_prim:invalidarg', ...
        'rv should be a positive integer with rv <= n.');
end
rv = int32(rv);

if nargin < 4 || isempty(msk)
    msk = [];
else
    if ~(islogical(msk) && isvector(msk) && numel(msk) == n)
        error('gr_prim:invalidarg', ...
            'msk should be a logical vector of length n.');
    end
end

if nargin < 5 || isempty(max_size)
    max_size = -1;
else
    if ~is_pos_scalar(max_size) 
        error('gr_prim:invalidarg', ...
            'max_size should be a positive scalar.');
    end
    if max_size >= n
        max_size = -1;
    end
end
max_size = int32(max_size);

if nargin < 6 || isempty(max_w)
    max_w = -1;
else
    if ~is_pos_scalar(max_w) 
        error('gr_prim:invalidarg', ...
            'max_size should be a positive scalar.');
    end
    if isinf(max_w)
        max_w = -1;
    end
end
max_w = cast(max_w, class(w));


%% main

edges = gr_prim_cimp(G, w, rv, msk, max_size, max_w);


%% Auxiliary function

function tf = is_pos_scalar(x)

tf = isnumeric(x) && isscalar(x) && isreal(x) && x == fix(x) && x > 0;



