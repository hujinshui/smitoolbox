function R = rand_label(K, n, op)
% Generates random labeling
%
%   R = rand_label(K, n);
%   R = rand_label(K, n, 'v');
%
%       Generates a 1 x n row vector, with each element being a
%       random integer in [1, K].
%
%   R = rand_label(K, n, 'b');
%
%       Generates a K x n matrix, with each column being a random 
%       indicator vector (a vector with a random entry set to one,
%       and others set to zeros).
%
%   R = rand_label(K, n, 'p');
%
%       Generates a K x n random matrix, such that each column has 
%       a unit sum.
%

%   Created by Dahua Lin, on Sep 4, 2011
%


%% verify input

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('rand_label:invalidarg', ...
        'K should be a positive integer scalar.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 0)
    error('rand_label:invalidarg', ...
        'n should be a non-negative integer scalar.');
end

if nargin < 3
    op = 'v';
else
    if ~(ischar(op) && isscalar(op))
        error('rand_label:invalidarg', ...
            'The 3rd argument should be a char.');
    end
    op = lower(op);
end

%% main

if op == 'v'
    R = randi(K, 1, n);
    
elseif op == 'b'
    L = randi(K, 1, n);
    R = l2mat(K, L);    
    
elseif op == 'p'
    R = rand(K, n);
    R = bsxfun(@times, R, 1 ./ sum(R, 1));
    
else
    error('rand_label:invalidarg', 'The 3rd argument is invalid.');
end




