function [dists, preds] = gr_dijkstra(G, s)
% Solve single-source shortest path problem using Dijkstra's algorithm
%
%   [dists, preds] = gr_dijkstra(G, s);
%       solves the single-source shortest-path problem with source s,
%       using Dijkstra's algorithm.       
%
%       Output arguments:
%       - dists:    the distances from source [1 x nv weight_type]
%       - preds:    the parents of each node [1 x nv int32]
%                   The parent of the source node is set to 0,
%                   and that of the non-visited node is set to -1.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 10, 2010
%       - Modified by Dahua Lin, on Mar 31, 2011
%

%% verify input

if ~(isa(G, 'gr_adjlist') && G.is_weighted)
    error('gr_dijkstra:invalidarg', 'G should be a (weighted) gr_adjlist object.');
end

if ~(isnumeric(s) && isscalar(s) && s >= 1 && s <= G.nv)
    error('gr_dijkstra:invalidarg', 's should be an integer scalar in [1, n]');
end

%% main

[dists, preds] = gr_dijkstra_cimp(G, int32(s-1));

preds(s) = 0;
un_visited = preds == (1:G.nv);
if any(un_visited)
    preds(un_visited) = -1;
    dists(un_visited) = inf;
end


