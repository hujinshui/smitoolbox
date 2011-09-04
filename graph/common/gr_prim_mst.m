function [preds, edges] = gr_prim_mst(G, r)
% Solve the minimum spanning tree (MST) using Prim's algorithm
%
%   preds = gr_prim_mst(G, r);
%   [preds, edges] = gr_prim_mst(G, r);
%
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
%       r is the root node of the MST, from which the solution starts.
%
%       In the output, 
%       - preds:    an 1 x nv int32 row vector, with preds(i) being the
%                   parent neighbor of vertex i. In particular, it has
%                   preds(r) == 0.
%
%       - edges:    an 1 x nv int32 row vector, with edges(i) being the
%                   index of the edge that connects vertex i and its
%                   parent. In particular, it has edges(r) == 0.
%

% Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

if ~(isa(G, 'gr_adjlist') && ~G.is_directed && G.is_weighted)
    error('gr_prim_mst:invalidarg', ...
        'G should be a (undirected & weighted) gr_adjlist object.');
end

if ~(isnumeric(r) && isscalar(r) && r >= 1 && r <= G.nv)
    error('gr_prim_mst:invalidarg', ...
        'r should be an integer scalar in [1, G.nv].');
end


%% main

ri = int32(r-1);

if nargout <= 1
    preds = gr_prim_cimp(G, ri);
else
    [preds, edges] = gr_prim_cimp(G, ri);
end

preds(r) = 0;

if nargout >= 2
    edges(r) = 0;
end

