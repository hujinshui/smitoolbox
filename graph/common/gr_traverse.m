function nds = gr_traverse(G, s, op)
% Traverse the nodes of a graph starting from a specified node
%
%   nds = gr_traverse(G, s, 'bf');
%   nds = gr_traverse(G, s, 'df');
%       traverse the nodes of a graph G starting from node s 
%       using either breadth-first order ('bf') or depth-first
%       order ('df').
%
%       G can be a weight matrix, or a mgraph struct, s is the
%       array of starting node indices. If s contains multiple
%       nodes, then it performs traverse starting from these
%       nodes sequentially. The nodes visited in previous 
%       traversal will not be visited again. 
%
%       nds is the a row vector of node indices in traversed order,
%       which contains all nodes that are connected at least one
%       node in s.
%       

% Created by Dahua Lin, on Oct 2, 2010
%

%% verify input

G = mgraph(G);

if ~isa(s, 'int32')
    s = int32(s);
end

if strcmp(op, 'bf')
    op = 'b';
elseif strcmp(op, 'df')
    op = 'd';
else
    error('gr_traverse:invalidarg', 'The 3rd argument is invalid.');
end

%% main

nds = gr_traverse_cimp(G, s, op);

