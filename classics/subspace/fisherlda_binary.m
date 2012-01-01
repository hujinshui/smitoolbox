function c = fisherlda_binary(X, L, varargin)
%FISHERLDA_BINARY Fisher's Linear Discriminant Analysis on two-classes
%
%   c = FISHERLDA_BINARY(X, L ...);
%       performs two-class linear discrminant analysis (LDA).
%
%       In the input, X is the data matrix, with each column being a
%       sample. Let the size of X be d x n, then there are n samples.
%       L should be a numeric or logical row vector of size 1 x n.
%       L(i) is the label for the sample X(:,i), which can be either
%       0 or 1.
%
%       The output c is a d x 1 discriminant coefficient vector,
%       which is orthogonal to the boundary plane.
%
%       One can specify further options when needed. Available 
%       options are listed below.
%
%       - 'reg':    the regularization coefficient, which would be
%                   added to the diagonal entries of the estimated
%                   common covariance. (default = 0)
%
%       - 'weights':  the weights of samples, which should be a
%                     row vector of size 1 x n. 
%                     (default = [], indicating all samples have the
%                      the same weight).
%

% Created by Dahua Lin, on Nov 22, 2010
%

%% verify input arguments

if ~(isfloat(X) && isreal(X) && ndims(X) == 2)
    error('fisherlda_binary:invalidarg', 'X should be a real matrix.');
end
n = size(X, 2);

if ~((islogical(L) || isnumeric(L)) && isequal(size(L), [1 n]))
    error('fisherlda_binary:invalidarg', ...
        'L should be a logical/numeric vector of size 1 x n.');
end
if ~islogical(L); L = logical(L); end

r = 0;
ws = [];

if ~isempty(varargin)
    onames = varargin(1:2:end);
    ovals = varargin(2:2:end);
    
    if ~(numel(onames) == numel(ovals) && iscellstr(onames))
        error('fisherlda_binary:invalidarg', ...
            'The option list is invalid.');
    end
    
    for i = 1 : numel(onames)
        
        name = onames{i};
        v = ovals{i};
        
        switch name
            case 'reg'
                if ~(isfloat(v) && isscalar(v) && v > 0)
                    error('fisherlda_binary:invalidarg', ...
                        'reg should be a real positive scalar.');
                end
                r = v;
            case 'weights'
                if ~(isfloat(v) && isequal(size(v), [1 n]))
                    error('fisherlda_binary:invalidarg', ...
                        'weights should be a 1 x n numeric vector.');
                end
                ws = v;
        end        
    end
end


%% main

% partition data

X0 = X(:, ~L);
X1 = X(:, L);

if isempty(X0)
    error('fisherlda_binary:invalidarg', 'No samples are with label 0.');
end

if isempty(X1)
    error('fisherlda_binary:invalidarg', 'No samples are with label 1.');
end

if isempty(ws)
    n0 = size(X0, 2);
    n1 = size(X1, 2);
else
    w0 = ws(~L);
    w1 = ws(L);
    t0 = sum(w0);
    t1 = sum(w1);
end

% estimate mean and covariance

if isempty(ws)
    mu0 = sum(X0, 2) * (1/n0);
    mu1 = sum(X1, 2) * (1/n1);
    
    C0 = X0 * X0' - (n0 * mu0) * mu0';
    C1 = X1 * X1' - (n1 * mu1) * mu1';
    Sigma = (C0 + C1) * (1 / (n0 + n1));
    
else
    mu0 = X0 * (w0' / t0);
    mu1 = X1 * (w1' / t1);
    
    C0 = X0 * bsxfun(@times, X0, w0)' - (t0 * mu0) * mu0';
    C1 = X1 * bsxfun(@times, X1, w1)' - (t1 * mu1) * mu1';
    Sigma = (C0 + C1) * (1 / (t0 + t1));
end

if r ~= 0
    Sigma = adddiag(Sigma, r);
end

% solve discriminant direction

c = Sigma \ (mu1 - mu0);
c = c / norm(c);


