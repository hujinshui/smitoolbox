function [M, L, info] = kmeans_std(X, M0, varargin)
%KMEANS_STD Standard K-means algorithm
%
%   [M, L] = KMEANS_STD(X, K, ...);
%   [M, L] = KMEANS_STD(X, M0, ...);
%   [M, L] = KMEANS_STD({X, w}, K, ...);
%   [M, L] = KMEANS_STD({X, w}, M0, ...);
%       
%       This function implements the standard K-means algorithm.
%
%       Input: 
%       - X:    The matrix of input samples. Each column in X corresponds
%               to one sample.
%
%       - w:    The weights of samples (when the samples are weighted)
%
%       - K:    The number of clusters. 
%
%       - M0:   The initial centers.
%   
%       Note that if K instead of M0 is input, the function will randomly
%       select K samples from X as the initial centers. One can also 
%       use the function kmseed to generate the initial centers. 
%
%       Output:
%       - M:    The resultant means, which is a matrix of size d x K,
%               where d is the space dimension, and K is the number of
%               clusters.
%
%       - L:    The assignment of samples to clusters. Suppose X contains
%               n samples, then L is a 1 x n row vector, whose values are
%               in {1, 2, ..., K}, and L(i) corresponds to X(:,i).
%
%
%       Options:
%       
%       The user can specify options (in name/value pairs) to 
%       control the algorithm. The options are listed as below:
%
%       - MaxIter:  the maximum number of iterations (default = 100)
%
%       - TolC:     the maximum allowable number of label changes at
%                   convergence. (default = 0)
%
%       - TolFun:   the maximum allowable change of objective value at
%                   convergence (default = 1e-8)
%
%       - Display:  the level of displaying (default = 'off')
%                   - 'off':    display nothing
%                   - 'final':  display a brief summary when the procedure
%                               finishes
%                   - 'iter':   display at each iteration
%                  
%       The user can also use kmopts function to construct
%       an option struct and use it as an input argument in the place
%       of the name/value pairs.
%
%
%   [M, L, info] = KMEANS_STD( ... );
%
%       additionally returns the information of the K-means procedure.
%
%       info is a struct with the following fields:
%       - niters:       The number of elapsed iterations
%       - objv:         The objective value at last step (total cost)
%       - converged:    Whether the procedure converged
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%       - Modified by Dahua Lin, on Mar 27 2012
%
    

%% parse and verify input

if isnumeric(X)     
    wx = [];
elseif iscell(X) && numel(X) == 2
    wx = X{2};
    X = X{1};
end

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('kmeans_std:invalidarg', 'X should be a real matrix.');
end
[d, n] = size(X);

if ~isempty(wx)
    if ~(isfloat(wx) && isreal(wx) && isvector(wx) && numel(wx) == n)
        error('kmeans_std:invalidarg', ...
            'w should be a real vector of length n.');
    end    
    if size(wx, 2) > 1
        wx = wx.';    % turn w into a column vector
    end
end

if isscalar(M0)
    K = M0;
    M0 = [];
else
    K = size(M0, 2);
    if ~(isfloat(M0) && ndims(M0) == 2 && size(M0, 1) == d)
        error('kmeans_std:invalidarg', 'M0 should be a d x K numeric matrix.');
    end
end

if ~(isnumeric(K) && (K >= 1 && K <= n) && K == fix(K))
    error('kmeans_std:invalidarg', 'The value of K is invalid.');
end

% get options

opts = kmopts(varargin{:});

tolfun = opts.tolfun;
tolc = opts.tolc;
maxiter = opts.maxiter;
displevel = opts.displevel;


%% main

% initialize centers

if isempty(M0)    
    seeds = randpick(n, K); 
    M0 = X(:, seeds);
end
M = M0;

% initialize assignments and objective

cpre = calc_costs_pre(X);
costs = calc_costs(M, X, cpre);
[min_costs, L] = min(costs, [], 1);

objv = calc_objv(min_costs, wx);

% main loop

t = 0;
converged = false;

if displevel >= 2
    print_iter_header();
end


while ~converged && t < maxiter
    
    t = t + 1;
    
    % identify affected clusters
    
    if t == 1
        aff_cs = 1 : K;
    else        
        aff_cs = find(intcount(K, [L(chs), L_pre(chs)]));
    end
                      
    % update centers        
        
    to_rp = false(1, K);
    
    gs = intgroup(K, L);
    for k = aff_cs
        gk = gs{k};
        if ~isempty(gk)
            M(:,k) = calc_mean(X, wx, gk);
        else
            to_rp(k) = 1;   % indicate to repick the k-th center
        end
    end    
    
    % repick new centers if needed
    
    if any(to_rp)
        rps = find(to_rp);
        M(:, rps) = X(:, randpick(n, numel(rps)));
    end
        
    % re-compute costs
    
    if numel(aff_cs) >= K / 2
        costs = calc_costs(M, X, cpre);
    else
        costs(aff_cs, :) = calc_costs(M(:, aff_cs), X, cpre);
    end
    
    % re-assign samples to centers
    L_pre = L;
    [min_costs, L] = min(costs, [], 1);
    
    % determine convergence
    
    objv_pre = objv;
    objv = calc_objv(min_costs, wx);
    
    v_ch = objv - objv_pre;
    
    chs = find(L ~= L_pre);
    converged = (abs(v_ch) <= tolfun || numel(chs) <= tolc);
    
    % print iteration info
    
    if displevel >= 2
        print_iter(t, numel(chs), aff_cs, objv, v_ch);
    end    

end

if displevel >= 1
    print_final(t, objv, converged);
end

if nargout >= 3
    info.niters = t;
    info.objv = objv;
    info.converged = converged;
end



%% Auxiliary functions

function v = calc_mean(X, w, si)

if isempty(w)
    cn = numel(si);
    cw = constmat(cn, 1, 1/cn);
else
    cw = w(si);
    cw = cw * (1 / sum(cw));
end

v = X(:, si) * cw;


function cpre = calc_costs_pre(X)

cpre.sX2 = sum(X.^2, 1);


function C = calc_costs(M, X, cpre)

C = (-2) * (M' * X);
C = bsxfun(@plus, C, sum(M .^ 2, 1).');
C = bsxfun(@plus, C, cpre.sX2);



function v = calc_objv(cs, w)

if isempty(w)
    v = sum(cs);
else
    v = cs * w;
end


%% Printing functions

function print_iter_header()

fprintf(' Iter    # ch.assign (aff.clus)     objv (change)  \n');
fprintf('------------------------------------------------------\n');


function print_iter(it, ch, afc, objv, vch)

fprintf('%5d        %7d (%5d)   %12.6g (%.4g)\n', it, ch, numel(afc), objv, vch);


function print_final(it, objv, converged)

fprintf('K-means: total_cost = %.6g, converged = %d [total # iters = %d]\n', ...
    objv, converged, it);




