function [s, t, w] = gr_kruskal_mst(G)
% Solve the minimum spanning tree (MST) using Kruskal's algorithm
%
%   [s, t, w] = gr_kruskal_mst(G);
%       solves the minimum spanning tree (forest) of the graph G using
%       Kruskal's algorithm. Here, G should be a symmetric/undirected graph.
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
    error('gr_kruskal_mst:invalidarg', 'G should be a weighted graph');
end


%% main

if nargout <= 2
    [s, t] = gr_kruskal_cimp(G);
else
    [s, t, w] = gr_kruskal_cimp(G);
end


