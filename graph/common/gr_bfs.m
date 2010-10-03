function vs = gr_bfs(G, seeds)
% Performs BFS traversal
%
%   vs = gr_bfs(G, seeds);
%       performs breadth-first-search traversal along the graph G.
%       G can be an affinity matrix or a mgraph struct.
%
%       seeds is an array comprised of the indices of seed nodes.
%
%       Outputs:
%       - vs:   the sequence of nodes in visiting order 
%

% Created by Dahua Lin, on Oct 3, 2010
%

%% verify input

G = mgraph(G);

if ~(isnumeric(seeds) && ~issparse(seeds))
    error('gr_bfs:invalidarg', 'seeds should be a non-sparse numeric array.');
end

if ~isa(seeds, 'int32')
    seeds = int32(seeds);
end

%% main

vs = gr_bfs_cimp(G, seeds-1);

