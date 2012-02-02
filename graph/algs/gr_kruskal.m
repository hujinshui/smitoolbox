function [edges, ccs] = gr_kruskal(G, w, K)
%GR_KRUSKAL Kruskal's Algorithm for minimum spanning tree/forest
%
%   edges = GR_KRUSKAL(G, w);
%   edges = GR_KRUSKAL(G, w, K);
%
%       Finds the minimum spanning tree/forest of the given graph G.
%
%       Input arguments:
%       - G:        The input graph, which should be an undirected graph
%                   struct.
%
%       - w:        The edge weights (a vector of length G.m).
%
%       - K:        The number of trees(clusters) in the forest. 
%                   If K is omitted, it is by default set to 1, meaning
%                   to find the minimum spanning tree.
%
%                   Note: generally, K is the minimum allowable number of
%                   trees. If G in itself have K' > K connected components, 
%                   then K' trees will actually be found.
%
%       Output arguments:
%       - edges:    The list of edge indices in the spanning tree/forest.
%
%                   Let m be the number of all edges in G, and e = edges(i)
%                   be the index of the i-th edge found in edges, whose
%                   value ranges from 1 to 2*m.
%
%                   Note the edges are in non-decreasing order of edge
%                   weights.
%
%   [edges, ccs] = GR_KRUSKAL( ... );
%       
%       Additionally returns a cell array of clusters. In particular, 
%       ccs{k} is the vector of vertex indices in the k-th cluster.
%

% Created by Dahua Lin, on Jan 26, 2012
%

%% verify input argument

if ~(is_gr(G) && isequal(G.dty, 'u'))
    error('gr_kruskal:invalidarg', ...
        'G should be an undirected graph struct.');
end

if ~(isfloat(w) && isreal(w) && ~issparse(w) && numel(w) == G.m)
    error('gr_kruskal:invalidarg', ...
        'w should be a real vector of length G.m.');
end
w = [w w];

if nargin < 3
    K = int32(1);
else
    if ~(isnumeric(K) && isscalar(K) && isreal(K) && K == fix(K) && ...
            K >= 1 && K < G.n)
        error('gr_kruskal:invalidarg', ...
            'K should be a positive integer scalar with K < G.n.');
    end
    
    K = int32(K);
end


%% main

if nargout <= 1
    edges = gr_kruskal_cimp(G, w, K);
else
    [edges, ccs] = gr_kruskal_cimp(G, w, K);
end








