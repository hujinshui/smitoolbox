function ccs = gr_conncomp(G)
% Get connected components of an undirected graph
%
%   ccs = gr_conncomp(G);
%       gets the connected component of G. To guarantee correctness,
%       the caller should ensure that G is symmetric, i.e. if i->j
%       then j->i.
%
%       In the output, ccs is an 1 x c cell array, where each cell is
%       a row vector comprised of the indices of nodes in a connected
%       component.
%
%   Remarks
%   -------
%       The internal implementation is based on BFS.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 2, 2010
%

%% verify input

G = mgraph(G);

%% main

ccs = gr_conncomp_cimp(G);

