function Y = mvsum(As, X, W, op)
% Compute the (weighted) sum of matrix products
%
%   Y = mvsum(As, X, W);
%   Y = mvsum(As, X, W, op);
%       compute the (weighted sum) of products between matrices and
%       vectors, respectively given by As and X, according to the 
%       formulas below:
%
%       When op = 'N':
%           y = \sum_{i=1}^n  w(i) * A_i * x_i,
%       When op = 'T':
%           y = \sum_{i=1}^n  w(i) * A_i' * x_i,
%
%       Here, A_i is given by A(:,:,i) and x_i is given by X(:,i).
%       The weights are given in each row of W.
%
%       If W is empty or a row vector, then it produces a single column
%       vector y. W can also be a matrix that contains m rows, in which
%       case, it produces m different y vectors, each corresponding to 
%       a row in W. In particular, Y(:,k) is computed based on the
%       weights given in W(k,:).
%
%       If the dimension of x-space is d and that of the y-space is q,
%       then when op is 'N', As should be a q x d x n array. If op is
%       'T', As(:,:,i) gives A_i', and in this case, As should be of
%       size d x q x n.
%

%   History
%   -------
%       - Created by Dahua Lin, on Apr 5, 2010
%       - Modified by Dahua Lin, on Apr 6, 2010
%           - support multiple groups of weights.
%       - Modified by Dahua Lin, on Apr 15, 2010
%           - change error handling.
% 

%% parse and verify input arguments

if ~(isfloat(As) && ndims(As) <= 3) 
    error('mvsum:invalidarg', ...
        'As should be a numeric array with ndims(As) <= 3.');
end

if ~(isfloat(X) && ndims(X) == 2)
    error('mvsum:invalidarg', ...
        'X should be a numeric matrix.');
end

if ~(ischar(op) && isscalar(op) && (op == 'N' || op == 'T'))
    error('mvsum:invalidarg', 'op should be either ''N'' or ''T''.');
end

if op == 'N'
    [q, d, n] = size(As);
    if ~(size(X,1) == d && size(As,3) == n)
        error('mvsum:invalidarg', 'The size of As is inconsistent with that of X.');
    end
else
    [d, q, n] = size(As);
    if ~(size(X,1) == d && size(As,3) == n)
        error('mvsum:invalidarg', 'The size of As is inconsistent with that of X.');
    end
end

if ~isempty(W)
    if ~(isfloat(W) && ndims(W) == 2 && size(W,2) == n)
        error('mvsum:invalidarg', 'W should be a numeric matrix with n columns.');
    end
end



%% main

if op == 'N'
    Ae = reshape(As, [q, d * n]);
    WX = compute_wx(X, W);
    
    Y = Ae * WX;
    
else % op == 'T'
    I = reshape(1:q*n, q, n).';
    Ae = reshape(As(:, I(:)), d * n, q);
    WX = compute_wx(X, W);
    
    Y = Ae' * WX;
end
        

%% sub-functions

function WX = compute_wx(X, W)

m = size(W, 1);

if isempty(W)
    WX = X(:);
elseif m == 1
    WX = bsxfun(@times, W, X);
    WX = WX(:);
else
    [d, n] = size(X);   
    We = reshape(W.', 1, n * m);
    We = We(ones(d, 1), :);
    We = reshape(We, d * n, m);    
    WX = bsxfun(@times, X(:), We);
end
    
