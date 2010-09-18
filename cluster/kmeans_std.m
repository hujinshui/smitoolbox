function [labels, centers, converged, terr] = kmeans_std(X, K, varargin)
%KMEANS_STD Implements standard K-means clustering 
%
%   [labels, centers] = kmeans_std(X, K);
%       performs K-means clustering for the input column vectors in X,
%       using default options.
%
%       X should be a d x n numeric matrix, with each column representing
%       a sample. Here, d and n respectively represent space
%       dimensionality, and the number of input samples. K is the number
%       of clusters. It is required that K <= n. 
%
%       In the output, labels will be a 1 x n row vector, in which
%       labels(i) indicates which center the i-th sample is finally
%       assigned to. centers will be a d x K vector, with each column
%       representing a center. centers(:, i) is the i-th center.
%
%       K-means is a classic method for clustering in vector spaces. 
%       The procedure comprises following steps:
%           - initialize the K centers by randomly selecting K samples,
%             which are called "seeds", from the input data, or by other
%             mechanism.
%           - iterate the following cycle until convergence or that the 
%             maximum iteration number is reached.
%               - assign each sample to the nearest center
%               - update the centers (means) based on the new assignment
%
%       In K-means algorithm, convergence means no label re-assignment
%       at the final iteration.
%
%       It is possible that some centers are not assigned any samples
%       during the iterations. This function supports different strategies
%       to handle the case. The default strategy is to randomly pick a 
%       sample to replace that center.
%
%   [labels, centers] = kmeans_std(X, K, name1, value1, name2, value2, ...);
%       performs K-means clustering with specified options.
%
%       The user can specify options in form of name-value list following
%       the argument K. The options that is not specified in the list is
%       set to the default value implicitly.
%
%       Here is a list of available options:
%           - 'init':   specify how to select the initial centers, whose
%                       value can be one of the following strings:
%                       - 'randseed':   randomly select K-samples from the 
%                                       the input data as initial centers
%                                       (this is default behavior)
%                       - 'maxdist':    It selects the first center that is 
%                                       cloest to the overall mean, and then                                        
%                                       selects other centers that maximize
%                                       the minimum distances to the
%                                       selected centers.
%                       The option value can also be a d x K matrix, which
%                       directly gives the K initial centers as the columns
%                       of the matrix.
%
%           - 'nilact': The action taken to handle the "nil centers", i.e.
%                       the centers that are not assigned any samples
%                       at some iteration. The option value can be one of
%                       the following strings:
%                       - 'repick':     randomly pick a sample from the
%                                       input data to replace that center.
%                                       (this is the default behavior)
%                       - 'mdpick':     pick a sample that is farthest from
%                                       other centers to replace that
%                                       center.
%                       - 'error':      raise an error with identifier
%                                       'kmeans_std:nilcenter'
%                       - 'discard':    simply discard that center, as a
%                                       consequence, the final number of 
%                                       clusters may be less than K.
%
%           - 'maxiter': The maximum number of iterations. The iterations
%                        will be stopped when the maxiter is reached, even
%                        when it has not converged. (default = 200)
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
%
%   [labels, centers, converged] = kmeans_std( ... );
%       addtionally returns a logical value indicating whether the
%       procedure has converged.
%
%   [labels, centers, converged, terr] = kmeans_std( ... );
%       additionally returns the total squared distances between samples to
%       their corresponding centers. 
%
%       If consider K-means as coding or vector quantization, this value
%       reflects the total coding / quantization error. K-means is actually
%       an greedy optimization procedure to minimize this error (locally).
%
%   Remarks
%       - This function is designed to be easy-to-use, efficient, and
%         versatile.distance
%
%       - The only requirement to run this function is the basic MATLAB at
%         version above 7.4 (2007a). The function is self-contained, it
%         does not rely on any other m-files, mex-files, or toolboxes.
%
%         Hence, you can simply copy this file to any place that your
%         MATLAB can find and then enjoy it!
%
%       - The function adopts several strategies to improve the efficiency:
%           - reformulate the computation of L2-distance matrix based on
%             matrix multiplication, and thus take advantage of MATLAB's
%             great performace in matrix multiplication to significantly
%             boost the computation speed. (refer to the subfunction
%             compute_dists)
%           - Observing that in many of the iterations, only a very small
%             set subset of labels is affected, thus only a small portion
%             of the centers need to be updated. The function figures out
%             which centers need to be updated and only re-compute those
%             centers.
%           - The program structure is carefully designed. Though it looks
%             long, only a part of codes would run in each particular
%             setting. The flexibility is offered with only slight overhead
%             caused by several if-judgement.
%
%         In many cases, this function runs 10 - 100 times faster than the
%         original kmeans function in MATLAB.
%
%   Example
%
%       % Suppose there are 1000 samples in a 50-dimensional space, then 
%       % we can use a matrix X of size 50 x 1000 to represent the samples.
%
%       % cluster the samples in X into 10 clusters using default way
%
%       [labels, centers] = kmeans_std(X, 10);
%
%       % cluster the samples in X into 10 clusters, 
%       % and guarantee that the number of iterations does not exceed 100
%
%       [labels, centers] = kmeans_std(X, 10, 'maxiter', 100);
%
%       % cluster the samples, discarding the nil centers
%   
%       [labels, centers] = kmeans_std(X, 10, 'nilact', 'discard');
%
%       % cluster the samples starting from the initial centers given by C0
%       % C0 should be a d x K matrix.
%
%       [labels, centers] = kmeans_std(X, 10, 'init', C0);
%
%       % cluster the samples and show the iteration procedure
%
%       [labels, centers] = kmeans_std(X, 10, 'display', 'iter');
%
%       % you also want to know whether the procedure converges
%
%       [labels, centers, converged] = kmeans_std(X, 10);
%
%       % you want to get all information of the procedure for analysis
%
%       [labels, centers, converged, terr, info] = kmeans_std(X, 10);
%
%       % you can simultaneously specify multiple options, such as
%
%       [labels, centers] = kmeans_std(X, 10, 'maxiter', 100, 'display', 'iter');
%

%   History
%       - Created by Dahua Lin, on Jun 3, 2009
%       - Modified by Dahua Lin, on Jun 4, 2009
%           - fix small bugs
%       - Modified by Dahua Lin, on Oct 15, 2009
%           - adjust interfaces
%


%% parse and verify input arguments

error(nargchk(2, inf, nargin));

assert(isfloat(X) && ndims(X) == 2, 'kmeans_std:invalidarg', ...
    'X should be a numeric matrix of floating-point value type.');

assert(isscalar(K) && K > 0 && K == fix(K), 'kmeans_std:invalidarg', ...
    'K should be a positive integer scalar.');

[d, n] = size(X);

assert(K <= n, 'kmeans_std:invalidarg', ...
    'K should not exceed the number of input samples.');

% deal with options

% set default options
    
init_func = @pick_randseed;
centers = [];
nilact = 'repick';
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
        'kmeans_std:invalidsyntax', ...
        'The options are not correctly specified.');
            
    % parse and verify options
    
    for i = 1 : nopts
        
        name = names{i};
        value = values{i};
        
        switch name
            case 'init'
                if isnumeric(value)
                    assert(isfloat(value) && isequal(size(value), [d K]), ...
                        'kmeans_std:invalidoption', ...
                        'The initial centers should be a d x K numeric matrix.');
                    
                    centers = value;
                    
                elseif ischar(value)
                    switch value
                        case 'randseed'
                            init_func = @pick_randseed;                            
                        case 'maxdist'
                            init_func = @pick_maxdist;
                        otherwise
                            error('kmeans_std:invalidoption', ...
                                'Unknown initialization method %s', value);
                    end
                    
                else
                    error('kmeans_std:invalidoption', ...
                        'The option value for init should be either a d x K matrix or a string.');
                end
                
                                
            case 'nilact'
                assert(ischar(value), 'kmeans_std:invalidoption', ...
                    'The option value for nilact should be a string.');
                
                assert( ...
                    strcmp(value, 'repick') || ...
                    strcmp(value, 'mdpick') || ...
                    strcmp(value, 'error') || ...
                    strcmp(value, 'discard'), ...
                    'kmeans_std:invalidoption', ...
                    'Unknown nil-center action %s', value);
                
                nilact = value;
                
            case 'maxiter'
                assert(isscalar(value) && value == fix(value) && value > 0, ...
                    'kmeans_std:invalidoption', ...
                    'The option value for maxiter should be a positive integer.');
                
                maxiter = value;
                
            case 'display'
                assert(ischar(value), 'kmeans_std:invalidoption', ...
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
                        error('kmeans_std:invalidoption', ...
                            'The option value for display is invalid.');
                end
                
            otherwise
                error('kmeans_std:invalidoption', ...
                    'The option name %s is not supported.', name);
        
        end  % switch option-name
        
    end % for each option-pair
    
end % if has options

% indicators (switches) for items to output or display

want_iter_terr = (nargout >= 4 || disp_level >= DISP_ITER);
want_final_terr = (nargout >= 3 || disp_level >= DISP_FINAL);


%% main procedure

% initialization

% initialize centers
if isempty(centers)  % centers are not explicitly given    
    centers = init_func(X, K);
end

nclusters = size(centers, 2);


% iterations

converged = false;
niters = 0;

while niters < maxiter
    
    niters = niters + 1;
    
    if niters > 1
        prev_labels = labels;     
    end
    
    
    % E-step: assign samples to current centers
    % -----------------------------------------
    
    % compute distances to centers (dists is a K x n matrix)
    dists = compute_dists(centers, X);
    
    % make assignment and compute terr
    [mindists, labels] = min(dists, [], 1);
    
    
    % Compute required information and check convergence
    % --------------------------------------------------
    
    if niters == 1        
        nchanges = n;
    else
        reassigned = find(prev_labels ~= labels);
        nchanges = numel(reassigned);
    end
    
    if want_iter_terr
        terr = sum(mindists);
    end    
    
    if disp_level >= DISP_ITER  % display the info about current iter
        fprintf('Iter %d:  #clusters = %d,  #changes = %d, terr = %g\n', ...
            niters, nclusters, nchanges, terr);
    end    
        
    if nchanges == 0
        converged = true;
        break;
    end        
    
    
    % M-step: update centers (re-compute means)
    % -----------------------------------------
            
    % figure out which centers are affected by re-assignment
    if niters == 1
        affected_labels = 1 : nclusters;
    else
        affected_labels = union(prev_labels(reassigned), labels(reassigned));
    end

    % re-compute affected centers only (and deal with nil-centers)
    [centers, labelmap] = update_centers(centers, nclusters, ...
        X, labels, affected_labels, nilact);
    
    nclusters = size(centers, 2);

    % re-map the labels if necessary
    if ~isempty(labelmap)
        labels = labelmap(labels);
    end          
    
end


% final operations

% compute final terr (if needed)

if want_final_terr && ~want_iter_terr
    terr = sum(mindists);
end

% display information (if needed)

if disp_level >= DISP_FINAL
    
    if converged
        fprintf('K-means completed (converged) with terr = %g\n', terr);
    else
        fprintf('K-means completed (NOT converged) with terr = %g\n', terr);
    end
    
elseif disp_level >= DISP_NOTIFY && ~converged    
    
    fprintf('K-means does not converged!\n'); 
    
end



%% core functions

function dists = compute_dists(X0, X)
% compute pairwise distances between columns in X0 and X

dists = (-2) * (X0' * X);
dists = bsxfun(@plus, dists, sum(X0 .* X0, 1)');
dists = bsxfun(@plus, dists, sum(X .* X, 1));


function [centers, labelmap] = update_centers(centers, K, X, labels, affected_labels, nilact)

% group samples based on labels

G = group_samples(labels, K);

% re-compute centers (record those without assigned samples)

is_nil = false(1, K);

for i = affected_labels
    inds = G{i};
    if ~isempty(inds)
        centers(:, i) = compute_center(X(:, inds));
    else
        is_nil(i) = true;
    end
end
    
% deal-with nil-centers (if there are)

if any(is_nil)
    
    nils = find(is_nil);
    
    switch nilact
        case 'repick'
            centers(:, nils) = pick_randseed(X, numel(nils));
            labelmap = [];
            
        case 'mdpick'
            centers(:, nils) = pick_maxdist(X, numel(nils), centers);
            labelmap = [];
            
        case 'error'
            error('kmeans_std:nilcenter', ...
                'Encounter nil centers during K-means iterations.');
            
        case 'discard'
            centers(:, nils) = [];
            labelmap = zeros(1, K);
            labelmap(~is_nil) = 1 : size(centers, 2);             
    end
        
else
    labelmap = [];
end
    


%% auxiliary functions

function v = compute_center(X)
% compute the center (mean vector) of column vectors in X

v = sum(X, 2) / size(X, 2);


function S = pick_randseed(X, K)
% randomly pick K seeds from X 

n = size(X, 2);
if K == 1
    S = X(:, max(ceil(rand() * n), 1));
elseif K < n
    S = X(:, randsample(n, K));    
else
    S = X;
end


function S = pick_maxdist(X, K, S0)
% pick samples from X to based on maximin dist criteria

[d, n] = size(X);
S = zeros(d, K);

if nargin == 2  || isempty(S0)
    % select new samples with maximin dist from selected samples
   
    % first sample (maximin dist from others in X)
    overall_mean = sum(X, 2) / n;    
    dists = compute_dists(overall_mean, X);
    [md, y] = min(dists); %#ok<ASGLU>
    S(:, 1) = X(:, y);
    
    % other sample (maximin dist from selected)
    for i = 2 : K
        dists = compute_dists(S(:, 1:i-1), X);
        dists = min(dists, [], 1);
        [md, y] = max(dists); %#ok<ASGLU>
        S(:, i) = X(:, y);
    end            
    
else   
    % select new samples with maximin dist from union of selected samples and S0
    
    for i = 1 : K
        dists = compute_dists([S0, S(:, 1:i-1)], X);
        dists = min(dists, [], 1);
        [md, y] = max(dists); %#ok<ASGLU>
        S(:, i) = X(:, y);
    end
    
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


          