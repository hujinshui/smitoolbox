function G = make_gr_bnd(nbs)
% Make a degree-bounded graph struct
%
%   The degree-bounded graph struct is a compact representation of
%   a directed graph whose (outgoing) degree is upper bounded by K.
%   
%   The struct contains the following fields:
%   - tag:      a string: 'gr-bnd'
%   - n:        the number of nodes
%   - K:        the upper bound of outgoing degrees
%   - nbs:      a K x n matrix, where nbs(:, i) lists the outgoing 
%               neighbors of the i-th node. If the degree of this node
%               is less than K, some entries in nbs(:,i) can be set to 
%               zeros to indicate no connection. (type: int32)
%
%   The actual number of edges equals nnz(nbs)
%
%   G = make_gr_bnd(nbs);
%
%       constructs a degree-bounded graph struct given the neighbor matrix.
%

% Created by Dahua Lin, on Nov 30, 2011
%

%% verify input argument

if ~(isnumeric(nbs) && ndims(nbs) == 2 && isreal(nbs) && ~issparse(nbs))
    error('make_gr_nbs:invalidarg', ...
        'nbs should be a non-sparse real numeric matrix.');
end

%% main

if ~isa(nbs, 'int32')
    nbs = int32(nbs);
end
[K, n] = size(nbs);

G.tag = 'gr-bnd';
G.n = int32(n);
G.K = int32(K);
G.nbs = nbs;

