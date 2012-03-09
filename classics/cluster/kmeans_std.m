function [L, M, info] = kmeans_std(X, M0, varargin)
% Standard K-means algorithm
%
%   [L, M] = kmeans_std(X, K, ...);
%   [L, M] = kmeans_std(X, M0, ...);
%   [L, M] = kmeans_std({X, w}, K, ...);
%   [L, M] = kmeans_std({X, w}, M0, ...);
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
%       - MaxIter:  the maximum number of iterations (default = 100)
%
%       - TolC:     the maximum allowable number of label changes at
%                   convergence. (default = 0)
%
%       - Display:  the level of information displaying (default = 'off')
%                   It can take either of the following values:
%                   - 'off':    display nothing
%                   - 'iter':   display information for each iteration
%                   - 'final':  display final result information
%
%       - UcWarn:   whether to raise a warning if not converged.
%                   (default = false).
%
%       - Init:     the method to initialize centers. It can take either
%                   of the following values (default = 'km++'):
%                   - 'random':  randomly pick K distinct samples as centers
%                   - 'km++':    randomly pick K samples using Kmean++
%                   - 'mcinit':  randomly select the first one, and then
%                                sequentially select furthest samples 
%                                from selected centers.                   
%
%       - OnNil:    the action taken to deal with "nil clusters", the 
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
%       - DistFunc:  the distance function (or its name)
%                     It can take either of the following values:
%                     - 'sqL2':  use squared L2 distance as cost, and
%                                component-wise arithmetic mean as center
%                     - 'L1':    use L1 distance as cost, and 
%                                component-wise median as center
%                       
%                     Or, it can be a user-defined function handle. 
%                     (refer to the remarks below for details).
%
%       The user can also use kmeans_ex_set function to construct
%       an option struct and input it to this function.   
%
%   Remarks
%   -------
%       Briefly, compared to conventional implementation, this function
%       has additional features, as follows:
%
%       - It supports weighted samples. 
%
%       - It supports different methods in initialization and in 
%         handling the case where there are some centers without
%         assigned samples.
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
%       - It uses various ways to increase the run-time efficiency, such
%         as vectorizing the computation of distances, and only updating
%         the distances / centers to affected clusters at each iteration.
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
    error('kmeans_std:invalidarg', 'X should be a numeric matrix.');
end
[d, n] = size(X);

if ~isempty(wx)
    if ~(isfloat(wx) && isequal(size(wx), [1 n]) && isreal(wx))
        error('kmeans_std:invalidarg', 'w should be a real vector of size 1 x n.');
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
        error('kmeans_std:invalidarg', 'M0 should be a d x K numeric matrix.');
    end
end

if ~(isnumeric(K) && (K >= 1 && K <= n) && K == fix(K))
    error('kmeans_std:invalidarg', 'The value of K is invalid.');
end

% get options

if isempty(varargin)
    opts = kmeans_std_set();
else
    if isstruct(varargin{1})
        opts = kmeans_std_set(varargin{:});
    else        
        opts = kmeans_std_set([], varargin{:});
    end
end
    
dfunc = opts.distfunc;
costfunc = @(x, y) dfunc(x, y, 'c');


%% main

% initialize centers

if isempty(M0)    
    seeds = opts.initfunc(X, K, [], costfunc);    
    M0 = X(:, seeds);
end
M = M0;

% initialize assignment and status

am = 1:K;
[L, ss] = reassign_std(X, M, [], dfunc);
G = intgroup([1, K], L);


% main loop

it = 0;
converged = false;
on_nil_op = opts.onnilop;

dispLevel = opts.displevel;
if dispLevel >= 2
    print_iter_header();    
    print_iter(it, ss.costs, wx, am, n);
end


while ~converged && it < opts.maxiter
    
    it = it + 1;
    
    % update affected mean
    
    [M, any_nil, to_rps] = update_mean(X, wx, M, dfunc, am, G, on_nil_op);
                
    if any_nil && on_nil_op > 0   % repick nil centers
        rpi = find(to_rps);
               
        if isempty(wx)
            rpw = ss.costs;
        else
            rpw = ss.costs .* wx;
        end
        
        picked_inds = opts.rpickfunc(X, numel(rpi), rpw, costfunc);        
        M(:, rpi) = X(:, picked_inds);
    end
    
    % update assignment and status
    
    G_pre = G;
    L_pre = L;
  
    [L, ss] = reassign_std(X, M, L, dfunc, ss, am);
    G = intgroup([1, K], L);
    
    % identify affected means
    
    am = false(1, K);
    for k = 1 : K
        am(k) = ~isequal(G_pre{k}, G{k});
    end
    am = find(am);        
    
    % determine convergence
    
    ch = nnz(L ~= L_pre);
    converged = (ch <= opts.tolc);
                
    % display iteration info
    if dispLevel >= 2    
        print_iter(it, ss.costs, wx, am, ch);
    end
end


% display final info
if dispLevel >= 1
    print_final(it, ss, wx, converged);
end

if ~converged && opts.ucwarn
    warning('kmeans_std:unconverged', ...
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

if ~isempty(am)
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
                error('kmeans_std:nilcenter', 'Nil center encountered.');
            end
        end
    end
end


function [L, ss] = reassign_std(X, M, L, dfunc, ss, am)
% Do cost updating and re-assignment in standard way

if isempty(L)
    ss.type = 's';
    ss.cmat = dfunc(M, X, 'c');
else
    K = size(M, 2);
    if ~isempty(am)
        if numel(am) < 0.5 * K
            ss.cmat(am, :) = dfunc(M(:, am), X, 'c');
        else
            % when many am, don't bother to take sub-matrices
            ss.cmat = dfunc(M, X, 'c');
        end
    end
end

[ss.costs, L] = min(ss.cmat, [], 1);


%% Auxiliary functions


%% Printing functions

function print_iter_header()

fprintf(' Iter    # af.m.  # ch.ass.     # t.cost \n');
fprintf('-------------------------------------------\n');


function print_iter(it, costs, wx, am, ch)

if isempty(wx)
    tcv = sum(costs);
else
    tcv = dot(costs, wx);
end
fprintf('%5d    %7d   %7d  %12.6g\n', it, numel(am), ch, tcv);


function print_final(it, ss, wx, converged)

if isempty(wx)
    tcv = sum(ss.costs);
else
    tcv = dot(ss.costs, wx);
end
fprintf('K-means final: total_cost = %.6g, converged = %d [total # iters = %d]\n', ...
    tcv, converged, it);




