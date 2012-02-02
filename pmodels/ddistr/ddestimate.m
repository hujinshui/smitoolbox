function P = ddestimate(V, w, c0)
%DDESTIMATE Estimates discrete distribution(s)
%
%   P = ddestimate(V);
%   P = ddestimate(V, w);
%
%       performs maximum likelihood estimation of discrete distribution(s)
%       given (weighted) observations.
%
%       Input arguments:
%       - V:    the observed histograms, a matrix of size K x n.
%               Here, K is the number of classes, and n is the number of
%               observed histograms.
%
%               V can also be in form of a cell array as {K, z}, where
%               z is a row vector of size 1 x n, and each element in z
%               is an integer between 1 and K.
%
%       - w:    the weights of the histograms, a row vector of size 1 x n
%               or a matrix of size m x n (m groups of weights).
%
%       Output arguments:
%       - P:    the estimated distribution(s).
%               if w is omitted or a row vector, then P is a column vector
%               of size K x 1, if there are multiple groups of weights,
%               namely m > 1, then P is of size K x m.
%
%   P = ddestimate(V, [], c0);
%   P = ddestimate(V, w, c0);
%
%       performs maximum a posteriori estimation of discrete distribution.
%       
%       Here, the prior is given by c0, the prior counts, which can be
%       either a scalar, or a vector of size K x 1. If Dirichlet prior is
%       used, then c0 here equals alpha - 1.
%

% Created by Dahua Lin, on Dec 25, 2011
%

%% parse input

if nargin < 2
    w = [];
end

if nargin >= 3 && ~isempty(c0)
    if ~(isfloat(c0) && isreal(c0) && ...
            (isscalar(c0) || isequal(size(c0), [K 1])))
        error('accum_counts:invalidarg', ...
            'c0 should be either a scalar or a real vector of size K x 1.');
    end
else
    c0 = 0;
end

%% main

H = accum_counts(V, w);

if ~isequal(c0, 0)
    if isscalar(c0) || size(H, 2) == 1
        H = H + c0;
    else
        H = bsxfun(@plus, H, c0);
    end
end

P = bsxfun(@times, H, 1 ./ sum(H, 1));
  
