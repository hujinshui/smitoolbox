function s = kmseeds(data, K, method, varargin)
%KMSEEDS Selecting seeds for clustering
%
%   s = KMSEEDS({X}, K, method, ...);
%   s = KMSEEDS({X, costfun}, K, method, ...);
%   s = KMSEEDS(D, K, method, ...);
%
%       Selects K seeds from a set of samples for K-means or K-medoid 
%       clustering. 
%
%       Input arguments:
%       - X:        The matrix comprising samples as columns
%
%       - costfun:  The function handle to evaluate costs between
%                   samples and centers.
%                   (If omitted, it is considered to be @pwsqL2dist).
%
%       - D:        The pairwise distance matrix [n x n]
%
%       - K:        The number of clusters.
%
%       - method:   The name of method used to select the seeds:
%
%                   - 'random':     Randomly select K seeds
%
%                   - 'km++':       K-means++ selection protocol
%                                   The probability being selected is
%                                   proportional to the square of min
%                                   distance to existing seeds
%
%                   - 'farthest':   Greedily select the sample with
%                                   greatest min-dist to existing seeds
%                                   as the next seed.
%
%       One can specify other options to control the selection, in form
%       of name/value pairs:
%       
%       - 'pre':        Pre-selected seed indices, a vector with K' indices.
%                       (default = [], indicating no pre-selected seeds)
%
%       - 'pre-vec':    Pre-selected seed vectors, a matrix composed of
%                       K' vectors that have been selected as seeds.
%
%       Note that is 'pre' is given, then 'pre-vec' is ignored. In
%       addition, 'pre-vec' can only be used when X is given.
%
%       Also, when there are pre-selected seeds, this function returns
%       K new seeds that are not in the pre-selected set.
%

% Created by Dahua Lin, on Mar 27, 2012
%

%% verify input arguments

if iscell(data) && (numel(data) == 1 || numel(data) == 2) 
    
    X = data{1};    
    if numel(data) == 1
        costfun = @pwsqL2dist;
    else
        costfun = data{2};
        if ~isa(costfun, 'function_handle')
            error('kmseeds:invalidarg', 'costfun must be a function handle.');
        end
    end

    if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
        error('kmseeds:invalidarg', 'X should be a real matrix.');
    end
    [d, n] = size(X);
    
    has_X = 1;
    
elseif isnumeric(data)
    
    D = data;
    
    if ~(isfloat(D) && isreal(D) && ndims(D) == 2 && size(D,1) == size(D,2))
        error('kmseeds:invalidarg', 'D should be a real squared matrix.');
    end
    n = size(D, 2);
    
    has_X = 0;
    
else
    error('kmseeds:invalidarg', 'The first argument is invalid.');
end

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1 && K <= n)
    error('kmseeds:invalidarg', 'K should be an integer in [1, n].');
end

if ~ischar(method)
    error('kmseeds:invalidarg', 'method should be a char string.');
end

pre = [];
pre_v = [];

if ~isempty(varargin)
    
    onames = varargin(1:2:end);
    ovals = varargin(2:2:end);
    
    for i = 1 : numel(onames)
        cn = onames{i};
        cv = ovals{i};
        
        switch lower(cn)
            case 'pre'
                if ~(isnumeric(cv) && (isvector(cv) || isempty(cv)))
                    error('kmseeds:invalidarg', ...
                        'pre should be a numeric vector.');
                end
                pre = cv;
                
            case 'pre-v'
                if ~has_X
                    error('kmseeds:invalidarg', ...
                        'pre-v only applies when X is given.');
                end
                
                if ~(isnumeric(cv) && isfloat(cv) && isreal(cv) && ...
                        ndims(cv) == 2 && size(cv, 1) == d)
                    error('kmseeds:invalidarg', ...
                        'The value given to pre-v is invalid.');
                end
                pre_v = cv;
        end        
    end
end





%% main

% preparation

if has_X    
    dfun = @(i) costfun(X(:,i), X);     
else
    dfun = @(i) D(i, :);
end

min_dists = [];

if ~isempty(pre)
    if has_X
        D0 = costfun(X(:, pre), X);
    else
        D0 = D(pre, :);
    end
    min_dists = min(D0, [], 1);
    
elseif ~isempty(pre_v)
    D0 = costfun(pre_v, X);
    min_dists = min(D0, [], 1);
end
    

% do selection

s = zeros(1, K);
    
switch method
    case 'random'
        if isempty(min_dists)
            s = randpick(n, K);
        else
            s0 = find(min_dists > 0);
            s = s0(randpick(numel(s), K));
        end                        
        
    case 'km++'
        s = do_progsel(dfun, n, K, min_dists, @kmpp_choose);
        
    case 'farthest'
        s = do_progsel(dfun, n, K, min_dists, @farthest_choose);
        
    otherwise
        error('kmseeds:invalidarg', 'Unknown method %s', method);
end


%% core functions

function s = do_progsel(dfun, n, K, min_dists, cfun)

s = zeros(1, K);

if isempty(min_dists)
    s(1) = randi(n, 1);
    min_dists = dfun(s(1));
else
    s(1) = cfun(min_dists);
    min_dists = min(min_dists, dfun(s(1)));
end

for k = 2 : K
    s(k) = cfun(min_dists);
    min_dists = min(min_dists, dfun(s(k)));
end

function i = kmpp_choose(dists)

p = dists * (1/sum(dists));
i = ddsample(p.', 1);
    
function i = farthest_choose(dists)

[~, i] = max(dists);
    
    
