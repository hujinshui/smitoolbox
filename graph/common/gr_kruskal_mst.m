function es = gr_kruskal_mst(G)
% Solve the minimum spanning tree (MST) using Kruskal's algorithm
%
%   es = gr_kruskal_mst(G);
%       solves the minimum spanning tree (forest) of the graph G using
%       Kruskal's algorithm. Here, G should be a symmetric/undirected graph.
%
%       The functions outputs a list of edge indices in row vector es.
%

% Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

G = gr_adjlist(G, 'u');

if isempty(G.w)
    error('gr_kruskal_mst:invalidarg', 'G should be a weighted graph');
end


%% main

es = gr_kruskal_cimp(G);

