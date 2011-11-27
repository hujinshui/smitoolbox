function P = cartprod(nums, op)
%CARTPROD Gets the cartesian product of number arrays
%
%   P = cartprod([n1, n2, ..., nK]);
%       computes the Cartesian product {1:n1} x {1:n2} x ... x {1:nK},
%       the results are returned by P, which is a K x (n1 x n2 x ... x nK)
%       array, with each column representing a combination of values.
%
%       This function can be used in cases involving combination 
%       traverse.
%
%       By default, the first element changes most slowly, and the
%       last element changes fastest, in the output.
%
%   P = cartprod([n1, n2, ..., nK], 'first-fast');
%       computes the Cartesian product. In the contrast to the 
%       default manner, using the 'first-fast' option, the first
%       element changes the fastest in the output.
%
%   Example
%       cartprod([2, 3, 4]) outputs
%           1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2
%           1 1 1 1 2 2 2 2 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3
%           1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4
%       
%       cardprod([2, 3, 4], 'first-fast') outputs
%           1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
%           1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 1 1 2 2 3 3 
%           1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4
%

%   History
%   -------
%       - Created by Dahua Lin, on May 30, 2008
%       - Modified by Dahua Lin, on Mar 22, 2010
%           - add the option for setting whether the first element changes
%             the fatest or the last element changes the fastest.
%           - refactor the program based on a new sub-function repvec             
%

%% parse and verify input

if isempty(nums)    
    P = [];
    return
else
    assert(isvector(nums) && isnumeric(nums) && ...
        all(nums > 0) && all(nums == fix(nums)), ...        
        'cartprod:invalidarg', ...
        'nums should be a vector of positive integers.');
    K = numel(nums);
end

if nargin >= 2
    assert(ischar(op) && strcmp(op, 'first-fast'), ...
        'cartprod:invalidarg', ...
        'The 2nd argument can only be ''first-fast'' if specified.');
    ff = true;
else
    ff = false;
end

        
% for K = 1,2,3, use specific simplified implementation
% for K > 3, use generic implementation

if K == 1
    P = 1 : nums(1);
    
elseif K == 2
    n1 = nums(1);
    n2 = nums(2);
    
    if ~ff
        P = [repvec(1:n1, n2, 1); repvec(1:n2, 1, n1)];
    else
        P = [repvec(1:n1, 1, n2); repvec(1:n2, n1, 1)];
    end
    
elseif K == 3
    n1 = nums(1);
    n2 = nums(2);
    n3 = nums(3);
    
    if ~ff
        P = [ ...
            repvec(1:n1, n2 * n3, 1); ...
            repvec(1:n2, n3, n1); ...
            repvec(1:n3, 1, n1 * n2)];
    else
        P = [ ...
            repvec(1:n1, 1, n2 * n3); ...
            repvec(1:n2, n1, n3); ...
            repvec(1:n3, n1 * n2, 1)];
    end
    
else % K > 3 (generic case)
    
    if size(nums, 1) > 1    % make it a row vector
        nums = nums.'; 
    end
    
    fcp = [1, cumprod(nums(1:K-1))];    % forward-cumprod:  [1, n1, n1 * n2, ...]
    rcp = [1, cumprod(nums(K:-1:2))];   % backward-cumprod: [1, n3, n3 * n2, ...]
    
    if ~ff
        irs = rcp(K:-1:1);      % irs: inner repeating times
        ors = fcp;              % ors: outer repeating times
    else
        irs = fcp;      
        ors = rcp(K:-1:1);
    end
    
    N = fcp(K) * nums(K);      % N: the total number
    P = zeros(K, N);
    for i = 1 : K
        P(i, :) = repvec(1:nums(i), irs(i), ors(i));
    end       
end


function r = repvec(v, m, n)
% v: a row vector of length k
% make a row vector of length m x k x n by repeating each element
% of v by m times, and then repeating the entire one by n times
%

if m > 1
    r = v(ones(m, 1), :);    
else
    r = v;
end

if n > 1
    r = r(:);
    r = r(:, ones(1, n));
end
   
if size(r, 1) > 1
    r = reshape(r, 1, numel(r));
end












    