function [labels, medoids, converged, tcost, info] = kmedoid_std(C, K, varargin)
%MLKMEDOID Implements standard K-medoid clustering
%
%   [labels, medoids] = kmedoid_std(C, K);
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
%       K is the number of clusters, which should be a positive integer
%       with K <= n.
%
%       The the output, labels will be a 1 x n row vector, in which
%       labels(i) indicates which cluster the i-th sample is assigned
%       to. medoids will be a 1 x K row vector, where medoids(i) is
%       the index of the sample that serves as the medoid of the i-th
%       cluster.
%
%       K-medoid is a classic method for clustering based on pairwise
%       relations. This function implements a batch version of the 
%       K-medoid algorithm, which can be briefly described as follows
%           - initialize the K medoids by randomly selecting K samples
%             or from other mechanism (input by user).
%           - iterate the following cycle until convergence or when the
%             maximum iteration number is reached
%               - assign each sample to the medoid with lowest cost
%               - update the medoids for each cluster, by selecting
%                 the one with lowest total cost for the cluster.
%
%       When no re-assignment happens during an iteration, the procedure
%       is considered to converge.
%
%   [labels, medoids] = kmedoid_std(X, K, name1, value1, name2, value2, ...);
%       performs K-medoid algorithm with specified options.
%
%       The user can specify options in form of name-value list following
%       the argument K. The options that is not specified in the list
%       is set to the default value implicitly.
%
%       Here is a list of available options:
%           - 'init':   specify how to select the initial medoids,
%                       whose value can be one of the following strings:
%                       - 'randseed':   randomly select K-samples from 
%                                       the input data as initial medoids.
%                       - 'maxdist':    selects the first medoid with 
%                                       lowest total cost for all samples,
%                                       and then gradually select the
%                                       medoids that maximize the min
%                                       cost to selected medoids.
%                       The option value can also be a 1 x K vector, which
%                       directly specifies the indices of the samples
%                       that serve as initial medoids.
%
%           - 'maxiter': The maximum number of iterations. The iterations
%                        will stop when the maxiter is reached even the
%                        procedure has not converged.
%
%           - 'display': What information is displayed in the procedure.
%                        The value can be any of the following strings:
%                        - 'off':       nothing is displayed (default)
%                        - 'iter':      displays information for each
%                                       iteration
%                        - 'final':     only displays the information of final
%                                       output
%                        - 'notify':    only displays a notification
%                                       message when the procedure does not
%                                       converge.
%
%   [labels, medoids, converged] = kmedoid_std( ... );
%       addtionally returns a logical value indicating whether the
%       procedure has converged.
%
%   [labels, medoids, converged, tcost] = kmedoid_std( ... );
%       additionally returns the total cost of representing the clusters
%       using the selected medoids. 
%
%   [labels, medoids, converged, terr, info] = kmedoid_std( ... );
%       additionally returns the iteration information through the 5th
%       output argument. info is a struct array, and info(i) is a struct
%       containing the information on the i-th iteration, which comprises
%       the following fields:
%           - 'nchanges':       the number of label re-assignment (changes)
%                               at the iteration.
%
%           - 'tcost':          the total cost of assigning each sample to 
%                               the corresponding clusters in the results.
%
%   Remarks
%       - The algorithm is for performing K-medoid in small and moderate
%         scale problems in which the entire cost matrix can be computed
%         and placed in the memory.
%
%   Example
%
%       % Suppose there are 1000 samples, then the cost matrix C
%       % should be a 1000 x 1000 matrix, that is pre-computed.
%
%       % cluster the samples into 10 clusters using default way
%
%       [labels, medoids] = kmedoid_std(X, 10);
%
%       % cluster the samples into 10 clusters, 
%       % and guarantee that the number of iterations does not exceed 100
%
%       [labels, medoids] = kmedoid_std(X, 10, 'maxiter', 100);
%
%       % cluster the samples starting from the initial medoids given by I0
%       % C0 should be a 1 x K row vector.
%
%       [labels, medoids] = kmedoid_std(X, 10, 'init', I0);
%
%       % cluster the samples and show the iteration procedure
%
%       [labels, medoids] = kmedoid_std(X, 10, 'display', 'iter');
%
%       % you also want to know whether the procedure converges
%
%       [labels, medoids, converged] = kmedoid_std(X, 10);
%
%       % you want to get all information of the procedure for analysis
%
%       [labels, medoids, converged, terr, info] = kmedoid_std(X, 10);
%
%       % you can simultaneously specify multiple options, such as
%
%       [labels, medoids] = kmedoid_std(X, 10, 'maxiter', 100, 'display', 'iter');
%

%   History
%       - Created by Dahua Lin, on Jun 4, 2008
%

%% parse and verify input arguments

error(nargchk(2, inf, nargin));

n = size(C,1);
assert(isnumeric(C) && ndims(C) == 2 && size(C,2) == n, ...
    'mlkmedoid:invalidarg', ...
    'The cost matrix C should be a numeric square matrix.');

assert(isscalar(K) && K > 0 && K == fix(K), 'mlkmedoid:invalidarg', ...
    'K should be a positive integer.');

assert(K <= n, 'mlkmedoid:invalidarg', ...
    'K should not exceed the number of input samples.');

% deal with options

init_func = @pick_randseed;
medoids = [];
maxiter = 200;

DISP_OFF = 0;
DISP_NOTIFY = 1;
DISP_FINAL = 2;
DISP_ITER = 3;

disp_level = DISP_OFF;

if ~isempty(varargin)
    
    % pre-process option list
    
    names = varargin(1:2:end);
    values = varargin(2:2:end);
    
    nopts = numel(names);
    assert(numel(values) == nopts && iscellstr(names), ...
        'mlkmedoid:invalidsyntax', ...
        'The options are not correctly specified.');
            
    % parse and verify options
    
    for i = 1 : nopts
        
        name = names{i};
        value = values{i};
        
        switch name
            case 'init'
                if isnumeric(value)
                    assert(isfloat(value) && isequal(size(value), [1 K]), ...
                        'mlkmedoid:invalidoption', ...
                        'The initial medoids should be a 1 x K vector.');
                    
                    medoids = value;
                    
                elseif ischar(value)
                    switch value
                        case 'randseed'
                            init_func = @pick_randseed;                            
                        case 'maxdist'
                            init_func = @pick_maxdist;
                        otherwise
                            error('mlkmedoid:invalidoption', ...
                                'Unknown initialization method %s', value);
                    end
                    
                else
                    error('mlkmedoid:invalidoption', ...
                        'The option value for init should be either a d x K matrix or a string.');
                end
                                
            case 'maxiter'
                assert(isscalar(value) && value == fix(value) && value > 0, ...
                    'mlkmedoid:invalidoption', ...
                    'The option value for maxiter should be a positive integer.');
                
                maxiter = value;
                
            case 'display'
                assert(ischar(value), 'mlkmedoid:invalidoption', ...
                    'The option value for display should be a string.');
                
                switch value
                    case 'off'
                        disp_level = DISP_OFF;
                    case 'iter'
                        disp_level = DISP_ITER;
                    case 'final'
                        disp_level = DISP_FINAL;
                    case 'notify'
                        disp_level = DISP_NOTIFY;
                    otherwise
                        error('mlkmedoid:invalidoption', ...
                            'The option value for display is invalid.');
                end
                
            otherwise
                error('mlkmedoid:invalidoption', ...
                    'The option name %s is not supported.', name);
        
        end  % switch option-name
        
    end % for each option-pair
    
end % if has options


% indicators (switches) for items to output or display

want_iter_tcost = (nargout >= 4 || disp_level >= DISP_ITER);
want_final_tcost = (nargout >= 3 || disp_level >= DISP_FINAL);

out_iter_info = (nargout >= 5);

if out_iter_info
    
    if isinf(maxiter)
        error('mlkmedoid:invalidoption', ...
            'maxiter should not be infinitely large when you want output info struct.');
    end
    
    % pre-allocate information structure
    info = repmat(struct('nchanges', [], 'tcost', []), maxiter, 1);
    
end



%% main

% initialization

% initialize medoids
if isempty(medoids)
    medoids = init_func(C, K);
end

nclusters = size(medoids, 2);

% iterations

converged = false;
niters = 0;

while niters < maxiter
    
    niters = niters + 1;
    
    if niters > 1
        prev_labels = labels;
    end
    
    % Re-assignment: assign samples to current medoids
    % ------------------------------------------------
    
    [mincosts, labels] = min(C(medoids, :), [], 1);
    
    % Compute required information and check convergence
    % --------------------------------------------------
    
    if niters == 1
        nchanges = n;
    else
        reassigned = find(prev_labels ~= labels);
        nchanges = numel(reassigned);
    end
    
    if want_iter_tcost
        tcost = sum(mincosts);
    end    
    
    if disp_level >= DISP_ITER  % display the info about current iter
        fprintf('Iter %d:  #changes = %d, tcost = %g\n', ...
            niters, nchanges, tcost);
    end    
    
    if out_iter_info
        info(niters).nchanges = nchanges;
        info(niters).tcost = tcost;
    end
        
    if nchanges == 0
        converged = true;
        break;
    end      
    
    
    % Model-update: update medoids for each cluster
    % ---------------------------------------------
    
    G = group_samples(labels, nclusters);
    
    for i = 1 : nclusters
        cidx = G{i};
        
        [mc, y] = min(sum(C(cidx, cidx), 2));
        medoids(i) = cidx(y); %#ok<AGROW>
    end
        
end

% final operations

% compute final terr (if needed)

if want_final_tcost && ~want_iter_tcost
    tcost = sum(mincosts);
end

% display information (if needed)

if disp_level >= DISP_FINAL
    
    if converged
        fprintf('K-medoids completed (converged) with tcost = %g\n', tcost);
    else
        fprintf('K-medoids completed (NOT converged) with tcost = %g\n', tcost);
    end
    
elseif disp_level >= DISP_NOTIFY && ~converged    
    
    fprintf('K-medoids does not converged!\n'); 
    
end

% prepare info output (if needed)

if out_iter_info
    
    if niters < length(info)
        info = info(1:niters);
    end
    
end

%% auxiliary functions

function S = pick_randseed(C, K)
% randomly pick K seeds from X 

n = size(C, 2);
if K == 1
    S = max(ceil(rand() * n), 1);
elseif K < n
    S = randsample(n, K)';    
else
    S = 1 : n;
end


function S = pick_maxdist(C, K)
% pick samples from X to based on maximin dist criteria

S = zeros(1, K);

% first sample (min cost to all)
[mc, y] = min(sum(C, 2));
S(1) = y;

% other sample (maximin cost from selected)
for i = 2 : K
    [mc, y] = max(min(C(S(1:i-1), :), [], 1));
    S(i) = y;
end
    

function G = group_samples(labels, K)
% group samples based on labels
% G{i} is the indices of the samples with label i

n = numel(labels);

[slabels, sidx] = sort(labels);
dp = find(diff(slabels));

sp = [1, 1 + dp];
ep = [dp, n];
ng = length(sp);

G = cell(1, K);

for i = 1 : ng
    l = slabels(sp(i));
    inds = sidx(sp(i) : ep(i));
    G{l} = inds;
end


