function [lens, vs, preds] = gr_dijkstra(G, w, s)
%GR_DIJKSTRA Single-source shortest paths 
%
%   [vs, lens] = GR_DIJKSTRA(G, w, s);
%   [vs, lens, preds] = GR_DIJKSTRA(G, w, s);
%
%       Find the shortest paths to all vertices of the graph from 
%       specific source using Dijkstra's algorithm.
%
%       Input arguments
%       - G:        The graph with neighborhood system. Both directed and
%                   undirected graphs are supported.
%
%       - w:        The edge weights (a vector of length m)
%       
%       - s:        The source, which can be either a single vertex, or
%                   multiple vertices (considered as a single source as
%                   a whole)
%       
%       Output arguments
%       - vs:       The list of all visited vertices. The vertices in vs 
%                   are in ascending order of shortest path lengths.
%
%       - lens:     The vector of corresponding shortest path lengths. 
%
%       - preds:    The vector of corresponding predecessors.
%
%       Specifically, lens(i) and preds(i) are respectively the shortest
%       path length and predecessor (along the shortest path) of the 
%       vertex i.
%
%   Remarks
%   -------
%       - The edge weights in w should all be non-negative for using 
%         Dijkstra's algorithm, otherwise incorrect results might be
%         returned.
%

% Created by Dahua Lin, on Jan 26, 2012
%

%% verify input arguments

if ~(is_gr(G) && G.has_nbs)
    error('gr_dijkstra:invalidarg', ...
        'G should be a graph with neighbor system.');
end

if ~(isfloat(w) && isreal(w) && ~issparse(w) && isvector(w) && length(w) == G.m)
    error('gr_dijkstra:invalidarg', ...
        'w should be a non-sparse real vector of length G.m.');
end

if G.dty == 'u'
    w = [w w];      % no matter whether w is column or row, this is ok
end
        
if ~(isnumeric(s) && isreal(s) && ~isempty(s) && ~issparse(s) && isvector(s))
    error('gr_dijkstra:invalidarg', ...
        's should be a non-empty numeric vector of indices.');
end
s = int32(s);


%% main

if nargout <= 1
    lens = gr_dijkstra_cimp(G, w, s);
elseif nargout == 2
    [lens, vs] = gr_dijkstra_cimp(G, w, s);
else
    [lens, vs, preds] = gr_dijkstra_cimp(G, w, s);
end
        


