function [I, J, W] = gr_kruskal_mst(G)
% Solve the minimum spanning tree (MST) using Kruskal's algorithm
%
%   [I, J, W] = gr_kruskal_mst(G);
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
%

% Created by Dahua Lin, on Oct 8, 2010
%

%% verify input

G = mgraph(G);

%% main

[I, J, W] = gr_mst_cimp(G, 'k');

