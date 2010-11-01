function [vs, parents, tf] = gr_bfs_trees(G, seeds)
% Construct a tree/forest by breadth-first-search on an undirected graph
%
%   [vs, parents, tf] = gr_bfs_trees(G, seeds);
%       constructs a BFS tree/forest that cover the connected components
%       involving all seeds. If seeds are omitted, it uses 1:n as seeds.
%
%       In the output, vs is the vertex traversal order, and parents
%       are the parents corresponding to vs.
%       tf indicates whether the input graph g is actually a
%       singly-connected graph.
%

% Created by Dahua Lin, on Oct 31, 2010.
%

%% verify input

G = gr_adjlist(G, 'u');

if nargin < 2
    seeds = 1 : G.n;
else
    if ~(isnumeric(seeds) && ~issparse(seeds))
        error('gr_bfs_trees:invalidarg', 'seeds should be a non-sparse numeric array.');
    end
end

if ~isa(seeds, 'int32')
    seeds = int32(seeds);
end

%% main

[vs, parents, tf] = gr_bfs_trees_cimp(G, seeds-1);


