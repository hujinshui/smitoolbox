function r = randpick(n, k, rstream)
% Random sample without replacement
%
%   r = randpick(n, k);
%   r = randpick(n, k, rstream);
%
%       randomly pick k distinct numbers from 1:n. In the output, r is
%       a vector of size 1 x k. (Note that it should have k <= n).
%       Note that r is sorted.
%
%       The user can specify the random number stream to use.  
%       

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%       - Modified by Dahua Lin, on Aug 10, 2011
%           - using an more efficient suffling algorithm
%

%% verify input

if ~(isscalar(n) && n == fix(n) && n > 0)
    error('randpick:invalidarg', 'n should be a positive scalar integer.');
end

if ~(isscalar(k) && k == fix(k) && k > 0 && k <= n)
    error('randpick:invalidarg', ...
        'k should be a positive scalar integer in [1, n].');
end

drs = 1; % use default random stream

if nargin >= 3 && ~isempty(rstream)
    if ~isa(rstream, 'RandStream')
        error('randpick:invalidarg', ...
            'rstream should be an object of class RandStream.');
    end
    drs = 0;
end

%% main

n = double(n);
k = double(k);

if k == 1    
    if drs
        r = randi(n);    
    else
        r = randi(rstream, n);
    end
    
elseif k == 2
    
    if drs
        r = randi(n, 2, 1);
    else
        r = randi(rstream, n, 2, 1);
    end
    
    while r(2) == r(1)
        if drs
            r(2) = randi(n);
        else
            r(2) = randi(rstream, n);
        end
    end    
    
else    
    if drs
        seed = uint32(randi(4294967295)); 
    else
        seed = uint32(randi(rstream, 4294967295));
    end
    
    if k * 100 < n
        r = randpick_cimp(n, k, 0, seed); % rejection sampling
    else
        r = randpick_cimp(n, k, 1, seed); % random shuffling
    end        
end
    

