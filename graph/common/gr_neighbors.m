function [nbs, ws] = gr_neighbors(G, op)
% Extract the neighbors of each node in a graph
%
%   nbs = gr_neighbors(G);
%   nbs = gr_neighbors(G, 'o');
%       extract the (outgoing) neighbors of all nodes. 
%
%       In the input, G is a weight matrix or an mgraph struct. 
%       In the output, nbs is a cell array, where nbs{i} is a row
%       vector of indices of the (outgoing) neighbors of i.
%       
%   nbs = gr_neighbors(G. 'i');
%       extract the incoming neighbors of all nodes.
%
%   [nbs, ws] = gr_neighbors( ... );
%       additionally returns the corresponding edge weights in form
%       of cell array.
%

%   History
%   -------
%       - Created by Dahua Lin, on Oct 1, 2010
%

%% verify input

G = mgraph(G);

if nargin < 2
    op = 'o';
else
    if ~(ischar(op) && numel(op) == 1)
        error('gr_neighbors:invalidarg', 'The 2nd argument should be a character.');
    end
end

%% main

if nargout <= 1
    nbs = gr_neighbors_cimp(G, op);
else
    [nbs, ws] = gr_neighbors_cimp(G, op);
end

n = numel(nbs);
for i = 1 : n
    nbs{i} = nbs{i} + 1;
end
    
