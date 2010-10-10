function [vs, preds, ford, dtime, ftime] = gr_dfs(G, seeds)
% Performs DFS traversal
%
%   vs = gr_dfs(G, seeds);
%   [vs, preds] = gr_dfs(G, seeds);
%   [vs, preds, ford] = gr_dfs(G, seeds);
%   [vs, preds, ford, dtime, ftime] = gr_dfs(G, seeds);
%
%       performs breadth-first-search traversal along the graph G.
%       G can be an affinity matrix or a mgraph struct.
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
%       - dtime:    the time-stamps of discovery
%       - ftime:    the time-stamps of finishing
%

% Created by Dahua Lin, on Oct 3, 2010
%

%% verify input

G = mgraph(G);

if ~(isnumeric(seeds) && ~issparse(seeds))
    error('gr_dfs:invalidarg', 'seeds should be a non-sparse numeric array.');
end

if ~isa(seeds, 'int32')
    seeds = int32(seeds);
end

%% main

nout = nargout;

if nout <= 1
    vs = gr_search_cimp(G, seeds-1, 'd');
elseif nout == 2
    [vs, preds] = gr_search_cimp(G, seeds-1, 'd');
elseif nout == 3
    [vs, preds, ford] = gr_search_cimp(G, seeds-1, 'd');
elseif nout <= 5
    [vs, preds, ford, dtime, ftime] = gr_search_cimp(G, seeds-1, 'd');
end

