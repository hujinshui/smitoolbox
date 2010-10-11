function [I, J, W] = gr_kruskal_mst(G)
% Solve the minimum spanning tree (MST) using Kruskal's algorithm
%
%   [I, J, W] = gr_kruskal_mst(G);
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
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
    [I, J] = gr_kruskal_cimp(G);
else
    [I, J, W] = gr_kruskal_cimp(G);
end


