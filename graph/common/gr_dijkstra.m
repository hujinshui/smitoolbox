function [vs, dists, preds] = gr_dijkstra(G, s)
% Solve single-source shortest path problem using Dijkstra's algorithm
%
%   [vs, dists, preds] = gr_dijkstra(G, s);
%       solves the single-source shortest-path problem with source s,
%       using Dijkstra's algorithm.       
%
%       Output arguments:
%       - vs:       the sequence of accessible nodes in the order of closing
%       - dists:    the corresponding distances from source
%       - preds:    the corresponding parents
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

if ~(isa(G, 'gr_adjlist') && G.is_weighted)
    error('gr_dijkstra:invalidarg', 'G should be a (weighted) gr_adjlist object.');
end

if ~(isnumeric(s) && isscalar(s) && s >= 1 && s <= G.nv)
    error('gr_dijkstra:invalidarg', 's should be an integer scalar in [1, n]');
end

%% main

[vs, preds, dists] = gr_dijkstra_cimp(G, int32(s-1));

