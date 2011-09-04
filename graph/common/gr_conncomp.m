function ccs = gr_conncomp(G)
% Get connected components of an undirected graph
%
%   L = gr_conncomp(G);
%       gets the connected component of G. To guarantee correctness,
%       the caller should ensure that G is symmetric, i.e. if i->j
%       then j->i.
%
%       The output is a label map of size 1 x n, which assigns each 
%       vertex a label indicating which component it belongs to.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 10, 2010
%

%% verify input

if ~(isa(G, 'gr_adjlist') && G.dtype == 'u')
    error('gr_conncomp:invalidarg', ...
        'G should be a (undirected) gr_adjlist object.');
end

%% main

ccs = gr_conncomp_cimp(G);

