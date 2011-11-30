function [G, w] = knng(D, K, varargin)
% Construct K-nearest-neighbor graph
%
%   [G, w] = knng(D, K, ...);
%
%       constructs a graph G by connecting each vertex to its K nearest
%       neighbors.
%
%       Input arguments:
%       - D:    The distances between vertices. D can be input in 
%               different forms:
%               - n x n dense matrix, D(i, j) represents the distance
%                 between the i-th and j-th vertices. 
%               - n x n sparse matrix, if D(i, j) is non-zero, then it
%                 represents the distance between the i-th and j-th
%                 vertices. 
%               - {g, w}: an graph struct together with edge weight/dist 
%                 values.
%
%       - K:    the (maximum) number of nearest neighbors for each vertex.
%
%       In addition, one can specify the following options in form of
%       name/value list.
%
%       - 's':      Can be 'min' or 'max', respectively indicating 
%                   whether minimum or maximum values implies nearest.
%                   default = 'min'.
%
%       - 'thres':  The threshold. Only those neighbors whose 
%                   corresponding values are below (for 'min') or
%                   above (for 'max') the threshold can be included.
%                   default = [], indicating no thresholding.
%
%       - 'sym':    Whether to make symmetric graph. default = false.
%                   If set to true, then when j is included as a neighbor
%                   of i, then i is added as a neighbor of j, even when
%                   it is not included based on the selection criteria.
%                   
%                   Note that if sym is set to true, the degree of each
%                   vertex is not necessarily bounded by K.
%
%       The output G is an object of class gr_adjlist. It is directed
%       when sym is false, and undirected when sym is true.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 13, 2010
%       - Modified by Dahua Lin, on Oct 28, 2011
%

%% verify input arguments

if isnumeric(D)
    if ~(ndims(D) == 2 && size(D,1) == size(D,2))
        error('knng:invalidarg', 'D should be a symmetric matrix.');
    end
    n = size(D, 1); 
    
    if issparse(D)
        [g, w] = gr_from_wmat(D, 'd', 'nbs');
        use_g = 1;
    else
        use_g = 0;
    end
    
elseif iscell(D) && numel(D) == 2
    
    g = D{1};
    w = D{2};
    
    if ~(g.dty == 'u' && g.has_nbs)
        error('knng:invalidarg', ...
            'The input graph should be undirected and have neighborhood.');
    end
    
    n = double(g.n);
    if ~(isnumeric(w) && isvector(w) && numel(w) == g.m)
        error('knng:invalidarg', ...
            'w should be a numeric vector of length m.');
    end
    
    use_g = 1;
else
    error('knng:invalidarg', 'The 1st argument is invalid.');
end

if ~(isnumeric(K) && isscalar(K) && K >= 1 && K <= n-1)
    error('knng:invalidarg', 'K should be a positive scalar in [1, n-1]');
end

opts = struct('s', 'min', 'thres', [], 'sym', false);
opts = parlist(opts, varargin{:});

if strcmp(opts.s, 'min')
    use_min = true;
elseif strcmp(opts.s, 'max')
    use_min = false;
else
    error('knng:invalidarg', 'Option s can only be ''min'' or ''max''.');
end

if ~(isempty(opts.thres) || (isnumeric(opts.thres) && isscalar(opts.thres)))
    error('knng:invalidarg', 'thres should be either empty or a numeric scalar.');
end

if ~(islogical(opts.sym) && isscalar(opts.sym))
    error('knng:invalidarg', 'sym should be a logical scalar.');
end

%% main

% extract neighboring edges

if use_g
    [s, t, w] = knng_cimp(g, w, double(K), use_min); 
else
    if strcmp(opts.s, 'min')
        rdv = inf;
    else
        rdv = -inf;
    end
    D(1:(n+1):n*n) = rdv;
    
    [w, t] = top_k(D, opts.s, K, 2);
    s = repmat((1:n).', 1, K); 
    
    s = s(:);
    t = t(:);
    w = w(:);
end    

    
% filter edges

ue = [];
if ~isempty(opts.thres)
    if use_min
        ue = find(w > opts.thres);
    else
        ue = find(w < opts.thres);
    end
end

if ~isempty(ue)
    s(ue) = [];
    t(ue) = [];
    w(ue) = [];
end


% make graph

if ~opts.sym
    G = make_gr('d', n, s, t, 'nbs');
else
    sa = [s; t];
    ta = [t; s];
    wa = [w; w];
    
    c = find(sa < ta);
    sa = sa(c);
    ta = ta(c);
    wa = wa(c);   
    
    [s, t, w] = make_unique(n, sa, ta, wa);    
    G = make_gr('u', n, s, t, 'nbs');
end
    

function [s, t, w] = make_unique(n, s, t, w)

e = sub2ind([n n], s, t);
[e, i] = unique(e); %#ok<ASGLU>

s = s(i);
t = t(i);
w = w(i);


    
        
    

