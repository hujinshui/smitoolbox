function [vs, parents, dists] = gr_bfs(G, seeds)
% Performs BFS traversal
%
%   vs = gr_bfs(G, seeds);
%   [vs, parents] = gr_bfs(G, seeds);
%   [vs, parents, dists] = gr_bfs(G, seeds);
%
%       performs breadth-first-search traversal along the graph G.
%       G can be an affinity matrix or a mgraph struct.
%
%       seeds is an array comprised of the indices of seed nodes, which
%       are all discovered initially.
%
%       Outputs:
%       - vs:       the sequence of nodes in visiting order 
%       - preds:    the parents of each node in the search tree
%       - dists:    the distances from the root of search tree
%

% Created by Dahua Lin, on Oct 3, 2010
%

%% verify input

G = gr_adjlist(G);

if ~(isnumeric(seeds) && ~issparse(seeds))
    error('gr_bfs:invalidarg', 'seeds should be a non-sparse numeric array.');
end

if ~isa(seeds, 'int32')
    seeds = int32(seeds);
end

%% main

if nargout <= 1
    vs = gr_bfs_cimp(G, seeds-1);
    
elseif nargout == 2
    [vs, parents] = gr_bfs_cimp(G, seeds-1);
    
elseif nargout == 3
    [vs, parents, dists] = gr_bfs_cimp(G, seeds-1);
end


