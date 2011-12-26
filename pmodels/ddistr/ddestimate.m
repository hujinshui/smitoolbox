function P = ddestimate(H, w, c0)
%DDESTIMATE Estimates discrete distribution(s)
%
%   P = ddestimate(H);
%   P = ddestimate(H, w);
%
%       performs maximum likelihood estimation of discrete distribution(s)
%       given (weighted) observations.
%
%       Input arguments:
%       - H:    the observed histograms, a matrix of size K x n.
%               Here, K is the number of classes, and n is the number of
%               observed histograms.
%
%               H can also be in form of a cell array as {K, z}, where
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
%   P = ddestimate(H, [], c0);
%   P = ddestimate(H, w, c0);
%
%       performs maximum a posteriori estimation of discrete distribution.
%       
%       Here, the prior is given by c0, the prior counts, which can be
%       either a scalar, or a vector of size K x 1. If Dirichlet prior is
%       used, then c0 here equals alpha - 1.
%

% Created by Dahua Lin, on Dec 25, 2011
%

%% verify input

if isnumeric(H)
    if ~(isfloat(H) && isreal(H) && ndims(H) == 2)
        error('ddestimate:invalidarg', 'H should be a real matrix.');
    end
    uhist = 1;
    [K, n] = size(H);
    
elseif iscell(H) && numel(H) == 2
    K = H{1};
    z = H{2};
    
    if ~(isscalar(K) && isreal(K) && K >= 1 && K == fix(K))
        error('ddestimate:invalidarg', 'K should be a positive integer.');
    end
    
    if ~(isnumeric(z) && isreal(z) && ndims(z) == 2 && size(z, 1) == 1)
        error('ddestimate:invalidarg', 'z should be a numeric row vector.');
    end
    
    uhist = 0;
    n = size(z, 2);    
end

if nargin >= 2 && ~isempty(w)
    if ~(isfloat(w) && isreal(w) && ndims(w) == 2 && size(w, 2) == n)
        error('ddestimate:invalidarg', ...
            'w should be a real matrix with size(w, 2) == n.');
    end
    has_w = 1;
else
    has_w = 0;
end

if nargin >= 3 && ~isempty(c0)
    if ~(isfloat(c0) && isreal(c0) && ...
            (isscalar(c0) || isequal(size(c0), [K 1])))
        error('ddestimate:invalidarg', ...
            'c0 should be either a scalar or a real vector of size K x 1.');
    end
else
    c0 = 0;
end


%% main

% summarize observations

if uhist
    
    if has_w
        s = H * w';
    else
        if n == 1
            s = H;
        else
            s = sum(H, 2);
        end
    end
    
else
    
    if has_w
        s = aggreg(w, K, z, 'sum');
    else
        s = intcount(K, z).';
    end
    
end

% incorporate prior

if ~isequal(c0, 0)
    if isscalar(c0) || size(s, 2) == 1
        s = s + c0;
    else
        s = bsxfun(@plus, s, c0);
    end
end

% normalize to get P

P = bsxfun(@times, s, 1 ./ sum(s, 1));

  
