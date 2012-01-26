function ccs = gr_conncomps(G)
%GR_CONNCOMPS Connected components of undirected graph
%
%   ccs = GR_CONNCOMPS(G);
%
%       Finds connected components of the graph G.
%
%       G should be an undirected graph struct (with neighbor system).
%
%       ccs is a 1 x K cell array (K is the number of components), and
%       ccs{k} is a row vector of indices of the vertices in the k-th
%       component.
%

% Created by Dahua Lin, on Jan 25, 2012
%

%% verify input arguments

if ~(is_gr(G) && isequal(G.dty, 'u') && G.has_nbs)
    error('gr_is_connected:invalidarg', ...
        'G should be an undirected graph struct (with neighbors).');
end

%% main

ccs = gr_conncomps_cimp(G);

