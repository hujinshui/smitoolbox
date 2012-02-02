function S = pca_std(X, pd, varargin)
%PCA_STD Standard principal component analysis
%
%   S = PCA_STD(X);
%   S = PCA_STD(X, pd, ...);
%
%       trains a standard principal component analysis model based on 
%       the data given as columns of X.
%
%       Suppose there are n samples in d-dimensional space, then X should 
%       be a d x n matrix with each column representing a sample.
%
%       The dimension of the principal subspace is determined in different 
%       ways, depending on the 2nd argument pd:
%
%       - if pd is a positive integer, then the dimension of the principal
%         subspace is pd. Note that in this case, pd should have 
%         pd <= min(d, n-1), where n is the number of columns in X.
%
%       - if 0 < pd < 1, then it determines a minimum subspace that 
%         preserves at least ratio pd of the total variance. 
%         (e.g if pd = 0.9, then it preserves 90% of the variance).
%
%       - pd can be [] (empty) or omitted, then it determines
%         a subspace that preserves 99.9% of the variance.
%         (eqivalent to setting pd to 0.999).
%
%       - pd can be a function handle that supports the
%         usage as follows:
%
%           d = pd(eigvals);
%
%         It takes as input a sorted (in descending order)
%         vector of eigenvalues, and returns the dimension
%         of the principal subspace.
%
%       In the output, S is a struct with the following fields:
%       - tag:          a string: 'pca-std'
%       - d:            the dimension of the input space
%       - dp:           the dimension of the principal subspace
%       - B:            the basis matrix [size: d x dp]
%       - pevs:         the eigenvalues of the principal subspace [1 x dp]
%       - res:          the residue energy (sum of residual eigenvalues)
%
%       In addition, one can specify the following options
%       in name/value pairs (optionally).
%
%       - 'method':     the method to train PCA model, which
%                       can be either of the following strings:
%                       - 'auto': automatically decide a proper
%                                 method. (default)
%                       - 'std':  the standard method, which
%                                 solves the eigen-problem of X*X'
%                       - 'svd':  the method based on SVD
%                                 decompostion.
%                       - 'trans': the method based on solves
%                                  the eigen-problem of X'*X,
%                                  which is efficient when n < d.
%
%       - 'weights':    The weights of samples. If this option
%                       is specified, then the model will be
%                       trained on weighted samples.
%
%                       The weights should be a 1 x n row
%                       vector, with weights(i) being the
%                       weight of the i-th sample.
%
%       - 'center':     The pre-computed center of the data
%                       in the input space.
%
%                       It can be either a d x 1 vector, or 0.
%
%                       If not specified, the sample mean will
%                       serve as the center.
%

%   History
%   -------
%       - Created by Dahua Lin, on May 30, 2008
%       - Modified by Dahua Lin, on Jun 6, 2008
%       - Modified by Dahua Lin, on Nov 20, 2010
%

%% parse and verify input arguments

% verify X
if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('pca_std:invalidarg', 'X should be a real matrix.');
end

[d, n] = size(X);

% verify pd

dim_ub = min(d, n-1);
if nargin < 2 || isempty(pd)
    pd = 0.999;
else
    if (isnumeric(pd) && isscalar(pd) && pd > 0) || ...
            isa(pd, 'function_handle')
        
        if isnumeric(pd) && pd >= 1
            if ~(pd == fix(pd) && pd <= dim_ub)
                error('pca_std:invalidarg', ...
                    'When pd >= 1, pd should be an integer and pd <= min(d, n-1)');
            end
        end
    else
        error('pca_std:invalidarg', 'pd is invalid.');
    end
end

[method, weights, cen] = parse_options(varargin);


%% pre-process samples

% shift-center
if isempty(cen)
    if isempty(weights)
        cen = sum(X, 2) / n;
    else
        cen = X * (weights' / sum(weights));
    end
end

if ~isequal(cen, 0)
    X = bsxfun(@minus, X, cen);
end

% weight samples
if isempty(weights)
    tw = n;
else
    tw = sum(weights);
    X = bsxfun(@times, X, weights);
end

%% perform spectral analysis

% decide method
if strcmp(method, 'auto')
    if n^2 * (d + n) < d^3
        method = 'trans';
    else
        method = 'std';
    end
end

% compute
switch method
    case 'std'
        M = X * X';
        M = 0.5 * (M + M');
        [U, D] = eig(M);
        clear M;
        
        devs = diag(D);
        clear D;
        
    case 'trans'
        M = X' * X;
        M = 0.5 * (M + M');
        [V, D] = eig(M);
        clear M;
        
        devs = diag(D);
        clear D;
        
        U = X * V;
        clear V;
        
        U = bsxfun(@times, U, 1./sqrt(sum(U.*U, 1)));
        
    case 'svd'
        [U, D] = svd(X, 'econ');
        devs = diag(D) .^ 2;
        clear D;
end

% re-arrange in descending order
[sevs, si] = sort(devs, 1, 'descend');
sevs = max(sevs(1:dim_ub), 0);
U = U(:, si(1:dim_ub));

sevs = sevs * (1 / tw);
tvar = sum(sevs);

% select principal subspace

if isnumeric(pd)
    if pd >= 1
        p = pd;
    else
        vr = cumsum(sevs) * (1 / tvar);
        p = find(vr >= pd, 1);
        if isempty(p); p = dim_ub; end
    end
else
    p = pd(sevs);
end

if p < dim_ub
    pevs = sevs(1:p);
    Up = U(:, 1:p);
    res = tvar - sum(pevs);
else
    pevs = sevs;
    Up = U;
    res = 0;
end

% make output

S.tag = 'pca-std';
[S.d, S.pd] = size(Up);
S.B = Up;
S.cen = cen;
S.pevs = pevs.';
S.res = res;



%% Sub functions

function [method, weights, cen] = parse_options(params)
% parse options

% default options

method = 'auto';
weights = [];
cen = [];

% parse options

if ~isempty(params)
    
    names = params(1:2:end);
    values = params(2:2:end);
    
    nopts = numel(names);
    if ~(numel(values) == nopts && iscellstr(names))
        error('pca_std:invalidsyntax', ...
            'The name-value list is invalid.');
    end
    
    for i = 1 : nopts
        
        name = lower(names{i});
        v = values{i};
        
        switch name
            case 'method'
                if ~(ischar(v) && ismember(v, {'auto','std','svd','trans'}))
                    error('pca_std:invalidoption', ...
                        'The method option value should be a string.');
                end
                method = v;
                
            case 'weights'
                if ~(isfloat(v) && isvector(v) && numel(v) == n)
                    error('pca_std:invalidoption', ...
                        'The weights should be a numeric vector of length n.');
                end
                if size(v, 1) > 1; v = v.'; end
                weights = v;
                
            case 'center'
                if ~(isfloat(v) && isreal(v) && ...
                        (isequal(v, 0) || isequal(size(v), [d 1])))
                    error('pca_std:invalidoption', ...
                        'The center should be either 0 or d x 1 real vector.');
                end
                cen = v;
                
            otherwise
                error('pca_std:invalidoption', ...
                    'Unknown option %s for estimating PCA', name);
        end
    end
end

            