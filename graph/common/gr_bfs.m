function [vs, preds, dists] = gr_bfs(G, seeds)
% Performs BFS traversal
%
%   vs = gr_bfs(G, seeds);
%   [vs, preds] = gr_bfs(G, seeds);
%   [vs, preds, dists] = gr_bfs(G, seeds);
%
%       performs breadth-first-search traversal along the graph G.
%       G can be an affinity matrix or a mgraph struct.
%
%       seeds is an array comprised of the indices of seed nodes, which
%       are all discovered initially.
%
%       Outputs:
%       - vs:       the sequence of nodes in visiting order 
%       - preds:    the predecessors of each nodes
%       - dists:    the distances (# edges) from the seeds
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

nout = nargout;

if nout <= 1
    vs = gr_search_cimp(G, seeds-1, 'b');
elseif nout == 2
    [vs, preds] = gr_search_cimp(G, seeds-1, 'b');
elseif nout == 3
    [vs, preds, dists] = gr_search_cimp(G, seeds-1, 'b');
end

