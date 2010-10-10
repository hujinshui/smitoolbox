function [I, J, W] = gr_prim_mst(G, r)
% Solve the minimum spanning tree (MST) using Prim's algorithm
%
%   [I, J, W] = gr_prim_mst(G);
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
%
%   [I, J, W] = gr_prim_mst(G, r);
%       solves the minimum spanning tree with root r. The resulting
%       edges will cover the connected component that contains r.
%

% Created by Dahua Lin, on Oct 8, 2010
%

%% verify input

if nargin < 2
    r = 0;
end
G = mgraph(G);

%% main

[I, J, W] = gr_mst_cimp(G, 'p', double(r));

