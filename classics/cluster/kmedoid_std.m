function [L, M, info] = kmedoid_std(C, M0, varargin)
%MLKMEDOID Implements standard K-medoid clustering
%
%   [L, M] = kmedoid_std(C, M0, ...);
%   [L, M] = kmedoid_std(C, K, ...);
%   [L, M, info] = kmedoid_std( ... );
%
%       performs K-medoid clustering based on the pairwise cost matrix C.
%
%       For a set of n samples, the cost matrix C should be an n x n
%       real matrix, with C(i, j) represents the cost of assigning j
%       to the clusters with medoid being the i-th sample. 
%
%       C can be pairwise distance matrix or any other costs that 
%       the user desire. Though C is symmetric in typical applications, 
%       this need not to be the case, asymmetric cost matrix is allowed.
%
%       M0 is the initial selection of medoids (in form of a vector 
%       of indices). Alternatively, the user can input the number of
%       medoids K, and let the function to select the initial medoids.
%       The option 'init' determines which method is used for initial
%       selection. By default, a method similar to Kmeans++ is used.
%
%       The the output, labels will be a 1 x n row vector, in which
%       labels(i) indicates which cluster the i-th sample is assigned
%       to. medoids will be a 1 x K row vector, where medoids(i) is
%       the index of the sample that serves as the medoid of the i-th
%       cluster.
%
%       K-medoid is a classic method for clustering based on pairwise
%       relations, which iterate between the following two steps 
%       until convergence.
%           - assign each sample to the medoid with lowest cost
%           - update the medoids for each cluster, by selecting
%             the one with lowest total cost for the cluster.    
%
%       In addition, one can specify the following options in form 
%       of name/value pairs (or a struct created by kmedoid_std_set)
%       to customize the algorithm.
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
%       - Rstream:  the random number stream to use, which can be either
%                   empty (using default) or a RandStream object.
%                   (default = []).
%
%   Remarks
%       - The algorithm is for performing K-medoid in small and moderate
%         scale problems in which the entire cost matrix can be computed
%         and placed in the memory.

%   History
%       - Created by Dahua Lin, on Jun 4, 2008
%       - Modified by Dahua Lin, on Sep 28, 2010
%

%% parse and verify input arguments

error(nargchk(2, inf, nargin));

if ~(isnumeric(C) && ndims(C) == 2 && size(C,1) == size(C,2))
    error('kmedoid_std:invalidarg', ...
        'The cost matrix C should be a numeric square matrix.');
end
n = size(C,1);

if isscalar(M0)
    K = M0;
    M0 = [];
    if ~(K == fix(K) && K >= 2 && K <= n) 
        error('kmedoid_std:invalidarg', 'K should be an integer in [2, n].');
    end    
else
    if ~(isvector(M0) && isnumeric(M0))
        error('kmedoid_std:invalidarg', 'M0 should be a vector of indices.');
    end
    if size(M0, 1) > 1
        M0 = M0.';
    end
    K = numel(M0);
end


% get options

if isempty(varargin)
    opts = kmedoid_std_set;
else
    if isstruct(varargin{1})
        opts = kmedoid_std_set(varargin{:});
    else        
        opts = kmedoid_std_set([], varargin{:});
    end
end


%% main

% initialization

% initialize medoids

if isempty(M0)
    M0 = opts.initfunc(n, K, [], C, opts.rstream);
end
M = M0;

% initialize assignment

cmap = C(M, :);
[costs, L] = min(cmap, [], 1);
G = intgroup([1, K], L);

% iterations

converged = false;
it = 0;
ac = 1 : K;

dispLevel = opts.displevel;
if dispLevel >= 2
    print_iter_header();    
    print_iter(it, costs, ac, n);
end

while ~converged && it < opts.maxiter
    
    it = it + 1;
    
    % update medoids
    
    for k = ac
        M(k) = select_medoid(C, G{k});
    end
    cmap(ac, :) = C(M(ac), :);
    
    % re-assign labels (only need to compare with the affected ones)
    
    G_pre = G;
    L_pre = L;
    
    [costs, L] = min(cmap, [], 1);    
    
    G = intgroup([1, K], L);
    
    % identify affected clusters
    ac = false(1, K);
    for k = 1 : K
        ac(k) = ~isequal(G_pre{k}, G{k});
    end
    ac = find(ac);
    
    % determine convergence
    
    ch = nnz(L ~= L_pre);
    converged = (ch <= opts.tolc);
                
    % display iteration info
    if dispLevel >= 2    
        print_iter(it, costs, ac, ch);
    end 
end


% display final info
if dispLevel >= 1
    print_final(it, costs, converged);
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



%% Sub functions


function i = select_medoid(C, g)

if isscalar(g)
    i = g;
else
    Cg = C(g, g);
    [mc, i] = min(sum(Cg, 2)); %#ok<ASGLU>
    i = g(i);
end


%% Printing functions

function print_iter_header()

fprintf(' Iter    # af.m.  # ch.ass.     # t.cost \n');
fprintf('-----------------------------------------------\n');


function print_iter(it, costs, am, ch)

tcv = sum(costs);
fprintf('%5d    %7d   %7d  %12.6g\n', it, numel(am), ch, tcv);


function print_final(it, costs, converged)

tcv = sum(costs);
fprintf('K-medoid final: total_cost = %.6g, converged = %d [total # iters = %d]\n', ...
    tcv, converged, it);




