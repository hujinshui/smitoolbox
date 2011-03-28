function G = ugprune(G0, K, w)
% Prune the edges of an undirected graph
%
%   G = ugprune(G0, K);
%   G = ugprune(G0, K, w);
%
%       makes a subgraph G from G0, by retaining up to K edges of largest 
%       weights for each vertex. 
%
%       If the 3rd argument is specified, the original weights are replaced
%       by the new weights in w. Here, w should be an 2m x 1 vector.
%
%       K can be empty. In this case, no edges will be removed.
%
%       Note that if an edge is within the K leading edges for one end
%       but not for the other, the edge is still retained.
%

% Created by Dahua Lin, on Feb 7, 2011
%

%% parse and verify input arguments

if ~isa(G0, 'gr_adjlist') && isequal(G0.dtype, 'u')
    error('ugprune:invalidarg', 'G0 should be an undirected gr_adjlist object.');
end

if ~(isempty(K) || (isscalar(K) && isnumeric(K) && K >= 1))
    error('ugprune:invalidarg', ...
        'K should be either empty or a positive integer number.');
end
K = double(K);

if nargin < 3
    w = G0.ew;
    if isempty(w)
        error('ugprune:invalidarg', 'The edge weights cannot be empty.');
    end
else
    m = G0.ne;
    if ~(isfloat(w) && isequal(size(w), [2*m, 1]))
        error('ugprune:invalidarg', ...
            'w should be a numeric vector of size 2m x 1.');
    end
end

if ~isa(w, 'double'); w = double(w); end


%% main

if isempty(K) || K >= max(G0.o_ds)  % no need to prune
    G = set_weights(G0, w);
    return;
end

% select edges

sel = ugprune_cimp(G0, K, w);
sel = sel + 1;

% make graph

s = G0.es(sel) + 1;
t = G0.et(sel) + 1;
w = w(sel);

G = gr_adjlist.from_edges('u', G0.nv, s, t, w);

