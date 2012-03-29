function r = randpick(n, k)
%RANDPICK Random sample without replacement
%
%   r = RANDPICK(n, k);
%
%       randomly pick k distinct numbers from 1:n. In the output, r is
%       a vector of size 1 x k. (Note that it should have k <= n).
%       Note that r is sorted.
%       

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%       - Modified by Dahua Lin, on Aug 10, 2011
%           - using an more efficient suffling algorithm
%       - Modified by Dahua Lin, on Jan 30, 2012
%

%% verify input

if ~(isscalar(n) && n == fix(n) && n > 0)
    error('randpick:invalidarg', 'n should be a positive scalar integer.');
end

if ~(isscalar(k) && k == fix(k) && k > 0 && k <= n)
    error('randpick:invalidarg', ...
        'k should be a positive scalar integer in [1, n].');
end

%% main

if k == 1    
    r = randi(n); 
    
else       
    n = int32(n);
    k = int32(k);
    
    use_shuffle = (n < k*200);
    r = randpick_cimp(n, k, use_shuffle);         
end
    

