function [G, V] = knng(D, op, K, thres)
% Construct K-nearest-neighbor graph
%
%   [G, V] = knng(D, 'min', K);
%   [G, V] = knng(D, 'max', K);
%   [G, V] = knng(D, 'min', K, thres);
%   [G, V] = knng(D, 'max', K, thres);
%       
%       constructs a K-nearest-neighbor graph given the
%       pairwise distances/similarities.
%
%       Inputs:
%       - D:        the n x n matrix of pairwise distances/similarities.
%                   Specifically, D(i, j) is the value between the node
%                   i and node j.
%
%                   The 2nd argument can be either 'min' or 'max'. If it
%                   is 'min', D represents distances (smaller value
%                   indicates closer); if it is 'max', D represents 
%                   similarities (greater value indicates closer).
%
%       - K:        the (maximum) number of neighbors of each node.
%
%       - thres:    If 2nd argument is 'min', only those edges whose
%                   associated distances are below thres are used.
%                   If 2nd argument is 'max', only those edges whose
%                   associated similarities are above thres are used.
%
%       Outputs:
%       - G:        the constructed graph, in form of a degree-bounded
%                   graph struct.
%       - V:        the corresponding edge distance/similarity matrix,
%                   of size K x n.
%
%       Note that the first row of G.nbs and w correspond to the nearest
%       neighbor, and the second row of G.nbs and w correspond to the
%       next nearest, and so on.
%
%       Moreover, loops (i.e. edges like (i, i)) are excluded.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 13, 2010
%       - Modified by Dahua Lin, on Oct 28, 2011
%       - Modified by Dahua Lin, on Nov 30, 2011
%

%% verify input arguments

if ~(isnumeric(D) && ndims(D) == 2 && ~issparse(D) && isreal(D) && ...
        size(D,1) == size(D,2))
    error('knng:invalidarg', 'D should be a non-sparse real square matrix.');
end
n = size(D,1);

if ischar(op)
    op = lower(op);
    if strcmp(op, 'max')
        is_max = true;
    elseif strcmp(op, 'min')
        is_max = false;
    else
        error('knng:invalidarg', ...
            'The 2nd argument should be either ''max'' or ''min''.');
    end    
else
    error('knng:invalidarg', ...
        'The 2nd argument should be either ''max'' or ''min''.');
end

if ~(isscalar(K) && isreal(K) && K == fix(K) && K >= 1 && K < n)
    error('knng:invalidarg', ...
        'K should be an integer in [1, n-1]');
end

if nargin < 4
    thres = [];
else
    if ~(isscalar(thres) && isreal(thres))
        error('knng:invalidarg', 'thres should be a real scalar.');
    end
end

%% main

Dr = rmdiag(D.', 'r');
[V, nbs] = top_k(Dr, op, K, 1);

a = bsxfun(@ge, nbs, 1:n);
nbs(a) = nbs(a) + 1;

if ~isempty(thres)
    if is_max
        nbs(V < thres) = 0;
    else
        nbs(V > thres) = 0;
    end
end

G = make_gr_bnd(nbs);

