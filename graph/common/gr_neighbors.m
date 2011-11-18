function C = gr_neighbors(G, vs)
% GR_NEIGHBORS Gets the neighbors of nodes
%
%   C = GR_NEIGHBORS(G);
%       Gets the neighbors of nodes of a graph and returns them as a cell
%       array.
%
%       Suppose G has n nodes, then C is a cell array of size n x 1, 
%       such that C{i} is a row vector of its neighbor node indices.
%
%   C = GR_NEIGHBORS(G, vs);
%       Gets the neighbors of the nodes whose indices are given in vs.
%

% Created by Dahua Lin, on Nov 10, 2011
%

%% verify input arguments

if ~is_gr(G)
    error('gr_neighbors:invalidarg', ...
        'G should be a struct with neighbors.');
end

if nargin < 2
    vs = 1:G.n;
else
    if ~(isnumeric(vs) && isreal(vs))
        error('gr_neighbors:invalidarg', 'vs should be a numeric array.');
    end
end

%% main

C = cell(size(vs));
os = G.o_os;
ds = G.o_degs;
nbs = G.o_nbs;

n = numel(vs);
for i = 1 : n    
    C{i} = nbs(os(i) + (1:ds(i)));        
end

