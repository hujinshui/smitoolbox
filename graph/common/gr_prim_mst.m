function es = gr_prim_mst(G, r)
% Solve the minimum spanning tree (MST) using Prim's algorithm
%
%   es = gr_prim_mst(G, r);
%
%       solves the minimum spanning tree (forest) of the graph G using
%       Prim's algorithm. Here, G should be a symmetric/undirected graph.
%       r is the root node of the MST, from which the solution starts.
%
%       The functions outputs a list of edge indices in row vector es.
%

% Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

if ~(isa(G, 'gr_adjlist') && ~G.is_directed && G.is_weighted)
    error('gr_prim_mst:invalidarg', ...
        'G should be a (undirected & weighted) gr_adjlist object.');
end


%% main

r = int32(r-1);
es = gr_prim_cimp(G, r);

m = G.ne;
rv = es > m;
es(rv) = es(rv) - m;


