function [vs, dists, preds] = gr_dijkstra_spath(G, s)
% Solve single-source shortest path problem using Dijkstra's algorithm
%
%   [vs, dists, preds] = gr_dijkstra_spath(G, s);
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
%       - Created by Dahua Lin, on Oct 4, 2010
%

%% verify input

G = mgraph(G);

if ~(isnumeric(s) && isscalar(s) && s >= 1 && s <= G.n)
    error('gr_dijkstra_spath:invalidarg', 's should be an integer scalar in [1, n]');
end

%% main

[vs, dists, preds] = gr_ss_spath_cimp(G, int32(s), 'd');

