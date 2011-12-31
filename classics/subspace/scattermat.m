function [mu, Sw, Sb] = scattermat(X, K, L, w)
% Compute scatter matrices for labeled Data
%   
%   [mu, Sw, Sb] = scattermat(X, K, L);
%   [mu, Sw, Sb] = scattermat(X, K, L, w);
%
%       computes within-class (and between-class) scatter matrices
%       from labeled data.
%
%       Inputs:
%       - X:    a d x n matrix, with each column being a sample.
%       - K:    the number of classes
%       - L:    the 1 x n label vector. L(i) is the label of sample X(:,i),
%               whose value should be in {1, 2, ..., K}.
%       - w:    the weights of samples (a vector of size 1 x n). 
%               If omitted, all samples have the same weight.
%
%       Outputs:
%       - mu:   the means of the samples
%       - Sw:   the within-class scatter matrix (i.e. pooled covariance)
%       - Sb:   the between-class scatter matrix (i.e. covariance of mean)
%

% Created by Dahua Lin, on Nov 22, 2010.
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('scattermat:invalidarg', 'X should be a real matrix.');
end

[d, n] = size(X);

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('scattermat:invalidarg', 'K should be a positive integer scalar.');
end

if ~(isnumeric(L) && isequal(size(L), [1 n]))
    error('scattermat:invalidarg', 'L should be a 1 x n numeric vector.');
end

if nargin < 4
    w = [];
else
    if ~(isempty(w) || (isfloat(w) && isequal(size(w), [1 n])))
        error('scattermat:invalidarg', 'w should be a 1 x n numeric vector.');
    end
end

%% main

% partition data

G = intgroup(K, L);

% compute mean

mu = zeros(d, K);
tw = zeros(1, K);

if isempty(w)
    for k = 1 : K
        gk = G{k};
        tw(k) = numel(gk);
        if tw(k) > 0
            mu(:,k) = sum(X(:,gk), 2) * (1/tw(k));
        end
    end
else
    for k = 1 : K
        gk = G{k};        
        if ~isempty(gk)
            wk = w(gk);
            tw(k) = sum(wk);
        else
            tw(k) = 0;
        end
        if tw(k) > 0
            mu(:,k) = X(:,gk) * (wk.' / tw(k));
        end
    end
end

wa = sum(tw);

% compute Sw

Z = X - mu(:, L);
if isempty(w)
    Sw = Z * Z';
else
    Sw = Z * bsxfun(@times, Z, w)';
    Sw = 0.5 * (Sw + Sw');
end

Sw = Sw * (1 / wa);

if nargin < 3
    return;
end

% compute Sb

cw = tw / wa;

mu0 = mu * cw';
Sb = mu * bsxfun(@times, mu, cw)' - mu0 * mu0';
Sb = 0.5 * (Sb + Sb');


