function [vs, preds, ford] = gr_dfs(G, seeds)
% Performs DFS traversal
%
%   vs = gr_dfs(G, seeds);
%   [vs, preds] = gr_dfs(G, seeds);
%   [vs, preds, ford] = gr_dfs(G, seeds);
%
%       performs depth-first-search traversal along the graph G.
%       G should be an object of class gr_adjlist.
%
%       seeds is an array comprised of the indices of seed nodes.
%       The search will be conducted as follows. It uses seeds(1)
%       for first DFS traversal, and then uses the next undiscovered
%       seed for next traversal.
%
%       Outputs:
%       - vs:       the sequence of nodes in visiting order 
%       - preds:    the predecessors of each nodes
%       - ford:     the finish order of nodes
%

% Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

if ~isa(G, 'gr_adjlist')
    error('gr_dfs:invalidarg', 'G should be a gr_adjlist object.');
end

if ~(isnumeric(seeds) && ~issparse(seeds))
    error('gr_dfs:invalidarg', 'seeds should be a non-sparse numeric array.');
end

if ~isa(seeds, 'int32')
    seeds = int32(seeds);
end

%% main

nout = nargout;

if nout <= 1
    vs = gr_dfs_cimp(G, seeds-1, 'd');
elseif nout == 2
    [vs, preds] = gr_dfs_cimp(G, seeds-1, 'd');
elseif nout == 3
    [vs, preds, ford] = gr_dfs_cimp(G, seeds-1, 'd');
end

