function W = cov_whiten(C, k, method, v)
% Solves a whitening transform
%
%   Given a Gaussian variable X with covariance C = U D U'.
%   Let W = D^{-1/2} * U', then WX has covarance I.
%   This transform W is called the whitening transform.
%   
%   W = cov_whiten(C);
%   W = cov_whiten(C, []);
%   W = cov_whiten(C, k);
%
%       computes the whitening transform W from covariance matrix C.
%       Note that the rows of W are sorted in descending order of norm.
%       
%       The second argument that only the first k rows are kepts in W.
%       If it is empty or omitted, then all rows are retained.
%
%   W = cov_whiten(C, k, 'bound', lb);
%   W = cov_whiten(C, k, 'reg', r);
%
%       In practice, it is often the case that C is singular or
%       ill-conditioned. To handle this case, we provide two approaches:
%
%       (1) bounding: bound the eigenvalues away from zero by force
%           all eigenvalues below thres = lb * max(eigenvalue)
%           to thres. If lb is omitted, default lb = 1e-3.
%
%       (2) regularization: add r * max(eigenvalue) to all diagonal
%           entries to C. If r is omitted, default r = 1e-3.
%
%   Remarks
%   -------
%       - The caller should ensure that the input matrix C is positive
%         semi definite.
%

% Created by Dahua Lin, on Nov 23, 2010
%

%% verify input argument

d = size(C, 1);
if ~(isfloat(C) && ndims(C) == 2 && isreal(C) && size(C,2) == d)
    error('cov_whiten:invalidarg', ...
        'C should be a real square matrix.');
end

if nargin < 2 || isempty(k)
    k = d;
else
    if ~(isnumeric(k) && isscalar(k) && k == fix(k) && k >= 1 && k <= d)
        error('cov_whiten:invalidarg', ...
            'k should be a positive integer scalar in [1, d].');
    end
end

if nargin < 3
    method = [];
else
    if ~ischar(method)
        error('cov_whiten:invalidarg', 'The 3rd argument should be a string.');
    end
    switch method
        case 'bound'
            if nargin < 4
                v = 1e-3;
            end
        case 'reg'
            if nargin < 4
                v = 1e-3;
            end
        otherwise
            error('cov_whiten:invalidarg', ...
                'The 3rd argument can only be ''bound'' or ''reg''.');
    end
end

%% main

[U, D] = eig(C);
evs = diag(D);

if ~isempty(method)
    switch method
        case 'bound'
            t = max(evs) * v;
            evs(evs < t) = t;
        case 'reg'
            r = max(evs) * v;
            evs = evs + r;
    end
end

w = 1 ./ sqrt(evs);
[w, si] = sort(w, 1, 'descend');

if k < d
    w = w(1:k);
    si = si(1:k);
end

U = U(:, si);
W = bsxfun(@times, w, U');


