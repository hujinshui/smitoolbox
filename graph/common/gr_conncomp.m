function ccs = gr_conncomp(G, op)
% Get connected components of an undirected graph
%
%   ccs = gr_conncomp(G);
%       gets the connected component of G. To guarantee correctness,
%       the caller should ensure that G is symmetric, i.e. if i->j
%       then j->i.
%
%       In the output, ccs is a cell array, where each cell is
%       a row vector comprised of the indices of nodes in a connected
%       component.
%
%   ccs = gr_conncomp(G, 'bf');
%   ccs = gr_conncomp(G, 'df');
%       specifies which traversal method to use. If omitted, breadth-first
%       search is used by default.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 2, 2010
%

%% verify input

G = mgraph(G);

if nargin < 2
    op = 'b';
else
    if strcmp(op, 'bf')
        op = 'b';
    elseif strcmp(op, 'df')
        op = 'd';
    else
        error('gr_traverse:invalidarg', 'The 3rd argument is invalid.');
    end
end

%% main

ccs = gr_conncomp_cimp(G, op);

