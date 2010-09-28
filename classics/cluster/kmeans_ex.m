function [L, M, info] = kmeans_ex(X, M0, varargin)
% Extended K-means algorithm
%
%   [L, M] = kmeans_ex(X, K, ...);
%   [L, M] = kmeans_ex(X, M0, ...);
%   [L, M] = kmeans_ex({X, w}, K, ...);
%   [L, M] = kmeans_ex({X, w}, M0, ...);
%       
%       This function implements the K-means algorithm, which is an
%       extension of the standard implementation.
%
%       Input: 
%       - X:    The matrix of input samples. Each column in X corresponds
%               to one sample.
%       - w:    The weights of samples (when the samples are weighted)
%       - K:    The number of clusters. 
%       - M0:   The initial centers.
%   
%       Note that if K instead of M0 is input, the function will invoke
%       a method to initialize the initial centers. The option 'init' 
%       controls how to do the initialization. By default, it uses km++ 
%       (Kmeans++) method.
%
%       Output:
%       - L:    The assignment of samples to clusters. Suppose X contains
%               n samples, then L is a 1 x n row vector, whose values are
%               in {1, 2, ..., K}, and L(i) corresponds to X(:,i).
%       - M:    The resultant means, which is a matrix of size d x K,
%               where d is the space dimension, and K is the number of
%               clusters.
%
%       Options:
%       
%       The user can specify options (in name/value pairs) to 
%       control the algorithm. The options are listed as below:
%
%       - scheme:   the scheme to use. (default = 'auto')
%                   It can be any of the following values:
%                   - 'std':   use standard scheme
%                   - 'acc':   use accelerated scheme
%
%       - max_iter: the maximum number of iterations (default = 100)
%
%       - tol_c:    the maximum allowable number of label changes at
%                   convergence. (default = 0)
%
%       - display:  the level of information displaying (default = 'off')
%                   It can take either of the following values:
%                   - 'off':    display nothing
%                   - 'iter':   display information for each iteration
%                   - 'final':  display final result information
%
%       - uc_warn:  whether to raise a warning if not converged.
%                   (default = false).
%
%       - init:     the method to initialize centers. It can take either
%                   of the following values (default = 'km++'):
%                   - 'random':  randomly pick K distinct samples as centers
%                   - 'km++':    randomly pick K samples using Kmean++
%                   - 'mcinit':  randomly select the first one, and then
%                                sequentially select furthest samples 
%                                from selected centers.                   
%
%       - on_nil:   the action taken to deal with "nil clusters", the 
%                   cluster to which no sample is assigned.
%                   (default = 'repick++').
%                   It can take either of the following values:
%                   - 'repick':     randomly repick a sample as new center
%                   - 'repick++':   randomly repick a sample as new center
%                                   using moderated scheme as in Kmeans++
%                   - 'mcpick':     pick the sample that is farthest to
%                                   all other centers (with maximum cost)
%                   - 'error':      raise an error.
%                   - 'keep':       simply keep that center without change.
%
%       - dist_func:  the distance function (or its name)
%                     It can take either of the following values:
%                     - 'sqL2':  use squared L2 distance as cost, and
%                                component-wise arithmetic mean as center
%                     - 'L1':    use L1 distance as cost, and 
%                                component-wise median as center
%                       
%                     Or, it can be a user-defined function handle. 
%                     (refer to the remarks below for details).
%
%       - rstream:  the random number stream to use, which can be either
%                   empty (using default) or a RandStream object.
%                   (default = []).
%
%       The user can also use kmeans_ex_set function to construct
%       an option struct and input it to this function.   
%
%   Remarks
%   -------
%       Briefly, compared to conventional implementation, this function
%       has additional features, as follows:
%
%       - It allows the use to choose between conventional implementation
%         and accelerated implementation. 
%
%         The conventional way computes the distances between all pairs of
%         samples and centers at every iteration, which leads to a great
%         deal of redundancy. The accelerated way reduces such redundant
%         computation exploiting triangle inequality, as described by the
%         following paper:
%
%           Charles Elkan. "Using the Triangle Inequality to Accelerate 
%           k-Means". Proceedings of 20th International Conference on
%           Machine Learning (ICML-03), Washington DC. 2003.
%
%         The basic idea of the accelerated method is to keep track of 
%         lower/upper bounds of relevant distances, and use triangle
%         inequality to identify unneccessary distance computation.
%         
%         We note that while the accelerated way avoids unnecessary 
%         computation, it would on the other hand incurs overhead of
%         book-keeping and subset selection. Moreover, the accelerated
%         scheme makes it difficult to vectorize some of the codes.
%         Hence, the overall run-time efficiency depends on particular
%         data set. By default, the function uses standard way.
%
%         Moreover, in our implementation, we avoid necessary
%         re-computation of mean-vectors, by only computing those affected
%         ones at each iteration.
%
%       - It allows the use of user-defined distance. The user can
%         implement its own distance by setting the dist_func option.
%         The function should support the following syntax:
%
%               R = dist_func(X, Y, 'd');
%                   computes the pairwise distances between X and Y.
%
%               R = dist_func(X, Y, 'c');
%                   computes the pairwise costs between X and Y. 
%                   Note that the cost used in objective function is not 
%                   necessary the same as distance. For example, in the
%                   standard Euclidean case, the cost of the square of
%                   the distance.
%
%               R = dist_func(X, [], 'm');
%               R = dist_func(X, w, 'm');
%                   computes the (weighted) mean of X.
%
%               R = dist_func(D, [], 't');
%                   computes the costs based on distances.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%
    

%% parse and verify input

if isnumeric(X)     
    wx = [];
elseif iscell(X) && numel(X) == 2
    wx = X{2};
    X = X{1};
end

if ~(isfloat(X) && ndims(X) == 2)
    error('kmeans_ex:invalidarg', 'X should be a numeric matrix.');
end
[d, n] = size(X);

if ~isempty(wx)
    if ~(isfloat(wx) && isequal(size(wx), [1 n]) && isreal(wx))
        error('kmeans_ex:invalidarg', 'w should be a real vector of size 1 x n.');
    end
    if issparse(wx)
        wx = full(wx);
    end
end

if isscalar(M0)
    K = M0;
    M0 = [];
else
    K = size(M0, 2);
    if ~(isfloat(M0) && ndims(M0) == 2 && size(M0, 1) == d)
        error('kmeans_ex:invalidarg', 'M0 should be a d x K numeric matrix.');
    end
end

if ~(isnumeric(K) && (K >= 1 && K <= n) && K == fix(K))
    error('kmeans_ex:invalidarg', 'The value of K is invalid.');
end

% get options

if isempty(varargin)
    opts = kmeans_ex_set;
else
    if isstruct(varargin{1})
        opts = kmeans_ex_set(varargin{:});
    else        
        opts = kmeans_ex_set([], varargin{:});
    end
end
    
dfunc = opts.dist_func;
costfunc = @(x, y) dfunc(x, y, 'c');

if strcmp(opts.scheme, 'std')
    use_acc = false;
else % 'acc'
    use_acc = true;
end


%% main

% initialize centers

if isempty(M0)    
    seeds = opts.initFunc(X, K, [], costfunc, opts.rstream);    
    M0 = X(:, seeds);
end
M = M0;

% initialize assignment and status

am = 1:K;
if use_acc
    [L, ss, dcc] = reassign_acc(X, M, [], [], [], dfunc);
else
    [L, ss, dcc] = reassign_std(X, M, [], dfunc);
end
G = intgroup([1, K], L);

tdcc = dcc;

% main loop

it = 0;
converged = false;
on_nil_op = opts.on_nil_op;

dispLevel = opts.dispLevel;
if dispLevel >= 2
    print_iter_header();    
    print_iter(it, ss.costs, wx, am, n, dcc);
end


while ~converged && it < opts.max_iter
    
    it = it + 1;
    
    % update affected mean
    
    M_pre = M;
    [M, any_nil, to_rps] = update_mean(X, wx, M, dfunc, am, G, on_nil_op);
                
    if any_nil && on_nil_op > 0   % repick nil centers
        rpi = find(to_rps);
               
        if isempty(wx)
            rpw = ss.costs;
        else
            rpw = ss.costs .* wx;
        end
        
        M(:, rpi) = opts.rpickFunc(X, numel(rpi), rpw, costfunc, opts.rstream);
    end
    
    % update assignment and status
    
    G_pre = G;
    L_pre = L;
    
    if use_acc
        [L, ss, dcc] = reassign_acc(X, M, M_pre, L, G, dfunc, ss, am);
    else
        [L, ss, dcc] = reassign_std(X, M, L, dfunc, ss, am);
    end    
    G = intgroup([1, K], L);
    tdcc = tdcc + dcc;
    
    % identify affected means
    
    am = false(1, K);
    for k = 1 : K
        am(k) = ~isequal(G_pre{k}, G{k});
    end
    am = find(am);
    
    % determine convergence
    
    ch = nnz(L ~= L_pre);
    converged = (ch <= opts.tol_c);
                
    % display iteration info
    if dispLevel >= 2    
        print_iter(it, ss.costs, wx, am, ch, dcc);
    end
end


% display final info
if dispLevel >= 1
    print_final(it, ss, wx, converged, tdcc);
end

if ~converged && opts.uc_warn
    warning('kmeans_ex:unconverged', ...
        'The K-means procedure did NOT converged.');
end
    
% output info
if nargout >= 3
    info.niters = it;
    info.last_ch = ch;
    info.converged = converged;
    
    if isempty(wx)
        info.totalcost = sum(ss.costs);
    else
        info.totalcost = dot(ss.costs, wx);
    end
end


%% core sub-functions


function [M, any_nil, to_rps] = update_mean(X, wx, M, dfunc, am, G, on_nil_op)
% update mean(s) according to new assignment
% and identify nil mean(s)

any_nil = 0;
if on_nil_op > 0
    to_rps = false(1, size(M, 2));
else
    to_rps = [];
end

for k = am    
    gk = G{k};
    if ~isempty(gk)        
        if isempty(wx)        
            M(:, k) = dfunc(X(:, gk), [], 'm');
        else
            M(:, k) = dfunc(X(:, gk), wx(gk), 'm');
        end
    else
        any_nil = 1;
        if on_nil_op > 0
            to_rps(k) = 1;
        elseif on_nil_op < 0
            error('kmeans_ex:nilcenter', 'Nil center encountered.');
        end
    end
end


function [L, ss, dcc] = reassign_std(X, M, L, dfunc, ss, am)
% Do cost updating and re-assignment in standard way

n = size(X, 2);

if isempty(L)
    ss.type = 's';
    ss.cmat = dfunc(M, X, 'c');
    dcc = n * size(M, 2);
else
    if ~isempty(am)
        ss.cmat(am, :) = dfunc(M(:, am), X, 'c');
        dcc = n * numel(am);
    else
        dcc = 0;
    end
end

[ss.costs, L] = min(ss.cmat, [], 1);


function [L, ss, dcc] = reassign_acc(X, M, M_pre, L, G, dfunc, ss, am)
% Do cost updating and re-assignment in accelerated way

n = size(X, 2);
K = size(M, 2);

if isempty(L)  % initialize
    
    ss.type = 'a';
    
    ss.DL = dfunc(M, X, 'd');   % lower-bound of pair-wise distances
    [ss.d, L] = min(ss.DL, [], 1);  % matching distances
    ss.costs = dfunc(ss.d, [], 't'); 
    
    ss.odated = false(size(ss.DL));  % true when DL(i,j) < dist(i, j)
        
    ss.Dm = dfunc(M, M, 'd');   % pairwise distances between means
    
    dcc = (n + K) * K;                
else
    
    if ~isempty(am)         
        
        % update center distances
        
        uDm = dfunc(M(:, am), M, 'd');
        ss.Dm(am, :) = uDm;
        ss.Dm(:, am) = uDm';
        
        % update lower bound
        
        for k = am
            mov_d = dfunc(M(:,k), M_pre(:,k), 'd');
            ss.DL(k, :) = max(ss.DL(k, :) - mov_d, 0);
        end        
        ss.odated(am, :) = true;    
                        
        dcc = numel(am) * (K + 1);
    else
        dcc = 0;
    end
    
    % update matching distances
    
    d = zeros(1, n);
    for k = 1 : K
        gk = G{k};
        if ~isempty(gk)
            d_gk = dfunc(M(:, k), X(:, gk), 'd');
            dcc = dcc + numel(gk);
            
            d(gk) = d_gk;
            ss.DL(k, gk) = d_gk;
            ss.odated(k, gk) = false;
        end
    end
    
    % update assignment (and thus new matching distances)
    
    for k = 1 : K
        
        % identify active candidate for re-assignment
        
        dl = ss.DL(k, :);
        hdc = ss.Dm(k, L) * 0.5;
        
        a = find(L ~= k & d > dl & d > hdc);

        if ~isempty(a)
                                
            da = ss.DL(k, a);
            
            % update active distances
            
            oa = ss.odated(k, a);
            oa_i = a(oa);
            
            if ~isempty(oa_i)                                                
                d_oa = dfunc(M(:, k), X(:, oa_i), 'd');
                dcc = dcc + numel(oa_i);
                
                da(oa) = d_oa;
                ss.DL(k, oa_i) = d_oa;
                ss.odated(k, oa_i) = false;
            end
            
            % do re-assignment
            j = find(da < d(a));            
            if ~isempty(j)
                da_j = da(j);
                j = a(j);
                
                L(j) = k;
                d(j) = da_j;
            end            
        end
                
    end
    
    ss.d = d;    
    ss.costs = dfunc(ss.d, [], 't');
end


%% Auxiliary functions


%% Printing functions

function print_iter_header()

fprintf(' Iter    # af.m.  # ch.ass.     # t.cost   # d.comp \n');
fprintf('-----------------------------------------------------\n');


function print_iter(it, costs, wx, am, ch, dcc)

if isempty(wx)
    tcv = sum(costs);
else
    tcv = dot(costs, wx);
end
fprintf('%5d    %7d   %7d  %12.6g  %9d\n', it, numel(am), ch, tcv, dcc);


function print_final(it, ss, wx, converged, tdcc)

if isempty(wx)
    tcv = sum(ss.costs);
else
    tcv = dot(ss.costs, wx);
end
fprintf('K-means final: total_cost = %.6g, converged = %d [total # iters = %d, # d.comp = %d]\n', ...
    tcv, converged, it, tdcc);




