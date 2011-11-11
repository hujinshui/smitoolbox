function [T, Tmu] = lda_ex(X, K, L, varargin)
% Extended Multi-class Linear Discriminant Analysis
%
%   T = lda_ex(X, K, L, ...);
%   [T, Tmu] = lda_ex(X, K, L, ...);
%       performs multi-class linear discriminant analysis. 
%
%       Suppose the training set comprises n samples on a d-dimensional
%       vector space. Then X should be a matrix of size d x n. Each
%       column of X is a sample.
%
%       K is the number of classes. L is a 1 x n row vector of the class 
%       labels. In particular, L(i) is the label for X(:,i), whose value 
%       should be in {1, 2, ..., K}, where K is the number of classes.
%
%       In the output, T is a transform matrix of size p x d (p < d), 
%       which can transform d-dimensional input vectors into 
%       p-dimensional discriminant features.
%
%       The 2nd output Tmu is a matrix of size p x K, where Tmu(:,k)
%       is the transformed center of the k-th class.
%
%       One can further specify some of the following options to 
%       control the algorithm in form of name/value pairs. 
%       For options that are not explicitly specified, the default
%       values will be used.
%
%       - 'maxdim': the maximum dimension of target space. 
%                   default = [], indicate to keep all dimensions.
%                   If this option is specified, at most maxdim dimensions
%                   are preserved.
%
%       - 'weights':  the sample weights. 
%                     default = [], indicate that all samples have the 
%                     same weight.
%
%       - 'wdim':   the dimension of whitened space. default = d.
%
%       - 'reg':    the regularization coefficient in whitening Sw.
%                   default = 1e-3 (i.e. add 1e-3 * max(eigenvalue)
%                   to diagonal entries of Sw before solving the
%                   whitening transform).
%
%       - 'bound':  the bounding coefficient in whitening Sw.
%                   default = [], indicate no bounding.
%                   If this option is specified to v, the function will
%                   enforce a lower bound = v * max(eigenvalue) to
%                   the eigenvalues of Sw before solving the whitening
%                   transform.
%
%       - 'ranktol':  the rank tolerance. The eigen-direction whose
%                     eigenvalue < ranktol * max(eigenvalue) is considered
%                     to be in the null space. default = 1e-8;
%
%   Remarks
%   -------
%       Here, we briefly introduce our implementation of LDA, in order
%       to show a clear picture of what the options mean, and how to
%       use them. The procedure has three main stages:
%
%       Stage 0: compute the mean vectors, and Scattering matrices,
%                including the within-class scatter matrix Sw, and
%                the between-class scatter matrix Sb. If the samples
%                are weighted, the weights will be utilized in the
%                computation in this stage.
%
%       Stage 1: Solve the whitening transform for Sw. In this stage,
%                either regularization or bounding can be done, depending
%                on whether reg or bound options are specified with
%                non-empty values. If they are simultaneously specified,
%                it applys bounding. In addition, we use wdim option
%                to control the dimension of whitened space.
%
%       Stage 2: After the whitening transform W is obtained, we compute
%                the matrix W * Sb * W', and solve the principal subspace
%                of the whitened space with respect to this matrix. 
%                The option ranktol (togther with maxdim) is used to 
%                determine the dimension of the transformed space.
%
%       Suppose the projection transform solved in stage 2 is P, then the
%       final transform that this function will return is P * W.
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2)
    error('lda_ex:invalidarg', 'X should be a numeric matrix.');
end
[d, n] = size(X);

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('lda_ex:invalidarg', 'K should be a positive integer scalar.');
end

if ~(isnumeric(L) && isequal(size(L), [1 n]))
    error('lda_ex:invalidarg', 'L should be a 1 x n numeric vector.');
end

% default options

maxdim = [];
w = [];
wdim = d;
reg = 1e-3;
bound = [];
rktol = 1e-8;

% parse options

if ~isempty(varargin)
    onames = varargin(1:2:end);
    ovals = varargin(2:2:end);
    if ~(numel(onames) == numel(ovals) && iscellstr(onames))
        error('lda_ex:invalidarg', 'The option list is invalid.');
    end
    
    for i = 1 : numel(onames)
        name = onames{i};
        v = ovals{i};
        switch lower(name)
            case 'maxdim'
                if ~(isnumeric(v) && isscalar(v) && v == fix(v) && v >= 1)
                    error('lda_ex:invalidopt', ...
                        'maxdim should be a positive integer scalar.');
                end
                maxdim = v;
            case 'weights'
                if ~(isempty(v) || (isfloat(v) && isequal(size(v), [1 n])))
                    error('lda_ex:invalidopt', ...
                        'weights should be either empty or an 1 x n numeric vector.');
                end
                w = v;
            case 'wdim'
                if ~(isnumeric(v) && isscalar(v) && v == fix(v) && v >= 1 && v <= d)
                    error('lda_ex:invalidopt', ...
                        'wdim should be an integer scalar in [1, d].');
                end
                wdim = v;
            case 'reg'
                if ~(isempty(v) || (isfloat(v) && isscalar(v) && isreal(v) && v > 0))
                    error('lda_ex:invalidopt', ...
                        'reg should be either empty or a positive scalar.');
                end
                reg = v;
            case 'bound'
                if ~(isempty(v) || (isfloat(v) && isscalar(v) && isreal(v) && v > 0))
                    error('lda_ex:invalidopt', ...
                        'bound should be either empty or a positive scalar.');
                end
                bound = v;
            case 'rktol'
                if ~(isempty(v) || (isfloat(v) && isscalar(v) && isreal(v) && v > 0))
                    error('lda_ex:invalidopt', ...
                        'ranktol should be either empty or a positive scalar.');
                end
                rktol = v;
            otherwise
                error('lda_ex:invalidarg', 'Unknown option name %s', name);
        end                
    end
end


%% main

% Stage 0: Compute Sw and Sb

[mu, Sw, Sb] = scattermat(X, K, L, w);

% Stage 1: Solve whitening transform w.r.t. Sw

if isempty(bound)
    if isempty(reg)
        W = cov_whiten(Sw, wdim);
    else
        W = cov_whiten(Sw, wdim, 'reg', reg);
    end
else
    W = cov_whiten(Sw, wdim, 'bound', bound);
end

% Stage 2: Solve projection w.r.t B = W * Sb * W'

B = W * Sb * W';
B = 0.5 * (B + B');

[U, D] = eig(B);
evs = diag(D);
[evs, si] = sort(evs, 1, 'descend');

if ~isempty(rktol)
    si = si(evs >= rktol * evs(1));
end

if numel(si) > maxdim;
    si = si(1:maxdim);
end
P = U(:, si)';

% Output

T = P * W;

if nargout >= 2
    Tmu = T * mu;
end


