function tf = gr_is_connected(G)
%GR_IS_CONNECTED Test graph connectivity
%
%   tf = GR_IS_CONNECTED(G);
%       Returns whether G is a connected graph. Here, G should be
%       an undirected graph.
%

% Created by Dahua Lin, on Jan 25, 2012
%

%% verify input arguments

if ~(is_gr(G) && isequal(G.dty, 'u') && G.has_nbs)
    error('gr_is_connected:invalidarg', ...
        'G should be an undirected graph struct (with neighbors).');
end

%% main

tf = gr_is_connected_cimp(G);

