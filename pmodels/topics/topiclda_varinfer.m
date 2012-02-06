function [Gamma, W, objVs, cvg] = topiclda_varinfer(Beta, alpha, C, w, Gamma, maxIter, tol)
%TOPICLDA_VARINFER Variational Inference Step of Latent Dirichlet Allocation
%
%   [Gamma, W] = TOPICLDA_VARINFER(Beta, alpha, C, w, Gamma, maxIter, tol);
%   [Gamma, W, objVs] = TOPICLDA_VARINFER( ... );
%   [Gamma, W, objVs, cvg] = TOPICLDA_VARINFER( ... );
%
%       Performs a variational inference update step for Latent Dirichlet
%       Allocation (LDA).
%
%       Suppose there are V words in vocabulary, K topics, and n documents.
%
%       Input arguments:
%       - Beta:     The word distributions of topics [V x K].
%
%       - alpha:    The Dirichlet prior parameter [scalar].
%
%       - C:        The word-count matrix [V x n]
%
%       - w:        The document-weights [empty or a vector of length n]
%
%       - Gamma:    The initial per-document topic distribution matrix [K x n]
%
%       - maxIter:  The maximum number of iterations (per document)
%
%       - tol:      The tolerance at convergence
%                   (in terms of L1-norm of gamma differnce)
%       
%       Output arguments:
%       - Gamma:    The solved per-document topic distribution matrix [K x n]
%
%       - W:        The accumulated per-topic word counts [K x V]
%
%       - objVs:    The per-document itemized objective term values [5 x n]
%                   Specifically, objVs(:,i) corresponds to the i-th doc.
%                   1st row:    expected log-likelihood of theta
%                   2nd row:    expected total log-likelihood of z
%                   3rd row:    expected total log-likelihood of w (words)
%                   4th row:    entropy of theta (w.r.t. gamma)
%                   5th row:    entropy of z (w.r.t. phi)
%
%       - cvg:      The vector of convergence indicator. [1 x n]
%                   cvg(i) indicates whether the inference on the i-th
%                   document converges.
%

% Created by Dahua Lin, on Feb 5, 2012
%

%% Verify inputs

if ~(isfloat(Beta) && isreal(Beta) && ndims(Beta) == 2 && ~issparse(Beta))
    error('topiclda_varinfer:invalidarg', ...
        'Beta should be a non-sparse real matrix.');
end
[V, K] = size(Beta);

if ~(isfloat(alpha) && isscalar(alpha) && isreal(alpha))
    error('topiclda_varinfer:invalidarg', 'alpha should be a real scalar.');
end
alpha = double(alpha);

if ~(isfloat(C) && isreal(C) && ndims(C) == 2 && size(C, 1) == V)
    error('topiclda_varinfer:invalidarg', ...
        'C should be a real matrix with V rows.');
end
n = size(C, 2);

if isempty(w)
    w = ones(1, n);
else
    if ~(isfloat(w) && isreal(w) && ~issparse(w) && isvector(w) && numel(w) == n)
        error('topiclda_varinfer:invalidarg', ...
            'w should be a real vector of length n.');
    end
    if ~isa(w, 'double'); w = double(w); end
end

if ~(isfloat(Gamma) && isreal(Gamma) && isequal(size(Gamma), [K n]))
    error('topiclda_varinfer:invalidarg', ...
        'Gamma should be a matrix of size K x n.');
end

if ~(isnumeric(maxIter) && isreal(maxIter) && isscalar(maxIter) && maxIter > 1)
    error('topiclda_varinfer:invalidarg', ...
        'maxIter should be a positive integer.');
end
maxIter = double(maxIter);

if ~(isnumeric(tol) && isreal(tol) && isscalar(tol) && tol > 0)
    error('topiclda_varinfer:invalidarg', ...
        'tol should be a positive real scalar.');
end
tol = double(tol);

    
%% main

[I, J, wc] = find(C);
I = int32(I) - 1;
J = int32(J) - 1; 
if ~isa(wc, 'double'); wc = double(wc); end

if ~isa(Gamma, 'double'); Gamma = double(Gamma); end

[Gamma, cvg, W, objVs] = ...
    topiclda_varinfer_cimp(Beta, alpha, w, I, J, wc, maxIter, tol, Gamma);



