function [L, cutv] = kolmogorov_mincut(g, caps, src_cap, snk_cap)
%KOLMOGOROV_MINCUT Kolmogorov's Mincut Algorithm
%
%   [L, cutv] = KOLMOGOROV_MINCUT(g, caps, src_cap, snk_cap);
%   [L, cutv] = KOLMOGOROV_MINCUT(g, caps, tcap);
%
%
%       Solve the mincut of a weighted graph, using Kolmogorov's mincut
%       algorithm.
%
%       Input arguments:
%       - g:                The graph stuct.
%
%                           If g is an undirected graph, then each edge
%                           is associated with capacities in both
%                           directions. 
%
%       - caps:             The vector of capacities of neighboring links
%                           (length = g.m)
%
%       - src_cap:          The capacities from source to all nodes
%                           (length = g.n)
%
%       - dst_cap:          The capacities from all nodes to sink 
%                           (length = g.n)
%
%       - tcap:             The signed capacities from/to terminals.
%                           (positive: from source, negative: to sink)
%
%       Output arguments:
%       - L:                The resulting label vector.
%                           L(i) = 1:    if node i is assigned to source
%                           L(i) = -1:   if node i is assigned to sink
%                           L(i) = 0:    node i can be assigned to either.
%
%       - cutv:             The total cut value of the minimum cut.
%

% Created by Dahua Lin, on Jan 29, 2011
%

%% verify input

if ~is_gr(g)
    error('kolmogorov_mincut:invalidarg', 'g should be a graph struct.');
end
n = double(g.n);
m = double(g.m);

if ~(isnumeric(caps) && isreal(caps) && isvector(caps) && numel(caps) == g.m)
    error('kolmogorov_mincut:invalidarg', ...
        'caps should be a real vector of length m.');
end

if nargin == 3
    tcap = src_cap;    
    
    if ~(isnumeric(tcap) && isreal(tcap) && isvector(tcap) && numel(tcap) == n)
        error('kolmogorov_mincut:invalidarg', ...
            'tcap should be a real vector of length n.');
    end
    
    src_cap(tcap < 0) = 0;
    
    snk_cap = zeros(size(tcap), class(tcap));
    snk_cap(tcap < 0) = -tcap(tcap < 0);        
else
    if ~(isnumeric(src_cap) && isreal(src_cap) && isvector(src_cap) && numel(src_cap) == n)
        error('kolmogorov_mincut:invalidarg', ...
            'src_cap should be a real vector of length n.');
    end    
    if ~(isnumeric(snk_cap) && isreal(snk_cap) && isvector(snk_cap) && numel(snk_cap) == n)
        error('kolmogorov_mincut:invalidarg', ...
            'snk_cap should be a real vector of length n.');
    end        
end

%% main

if g.dty == 'u'
    nb_cap = caps;
    rv_cap = caps;
else
    nb_cap = caps;
    rv_cap = zeros(size(caps), class(caps));
end

sv = g.edges(1, 1:g.m);
tv = g.edges(2, 1:g.m);

[L, cutv] = kolmogorov_mincut_cimp(n, m, sv, tv, nb_cap, rv_cap, src_cap, snk_cap);

