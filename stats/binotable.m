function M = binotable(n, p)
% Generate a table of binomial probabilities
%
%   M = binotable(n, p);
%       generates a table of size (n+1) x (n+1), such that
%       for each i = 0, ..., n, and each j = 0, ..., i,
%       M(i, j) is the probability of having j successes from i
%       independent Bernoulli trials with success rate p.
%
%       The formula is given by
%
%           M(i+1, j+1) = C(i, j) * p^j * (1 - p)^(i-j),
%
%       where, C(i, j) is the number of combinations with j items chosen
%       from i ones.
%
%       In the output, M(i, j) always equals zero when j > i. In other 
%       words, M is a lower triangle matrix.
%

% Created by Dahua Lin, on Apr 24, 2010
%

%% verify input arguments

if ~(isfloat(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('binotable:invalidarg', ...
        'n should be a non-negative integer scalar.');
end

if ~(isfloat(p) && isscalar(p) && isreal(p) && p >= 0 && p <= 1)
    error('binotable:invalidarg', ...
        'p should be a real scalar in [0, 1]');
end

%% main

if n == 0
    M = 1;
    
else % n > 0
    
    fv = [1 cumprod(1:n)];
    a = (p .^ (0:n)) ./ fv;
    b = ((1-p).^(0:n)) ./ fv;
    
    % get all triangle indices
    
    I = (0:n)' * ones(1, n+1);
    J = ones(n+1, 1) * (0:n);
    
    i = I(I >= J);
    j = J(I >= J);
    
    inds = (i+1) + j * (n+1);
    
    % generate the table
    
    M = zeros(n+1, n+1);
    M(inds) = fv(i+1) .* a(j+1) .* b(i-j+1);    
end
    
