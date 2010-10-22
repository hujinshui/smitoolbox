function [s, t, w] = gr_prim_mst(G, r)
% Solve the minimum spanning tree (MST) using Prim's algorithm
%
%   [s,t,w] = gr_prim_mst(G, r);
%
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
%       r is the root node of the MST, from which the solution starts.
%
%       In the output, 
%           - s:    the source vertices of edges
%           - t:    the target vertices of edges
%           - w:    the weights of edges
%

% Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

G = gr_adjlist(G);

if isempty(G.w)
    error('gr_prim_mst:invalidarg', 'G should be a weighted graph');
end


%% main

r = int32(r-1);

if nargout <= 2
    [s, t] = gr_prim_cimp(G, r);
elseif nargout == 3
    [s, t, w] = gr_prim_cimp(G, r);
end


