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
%   r = randpick(p, k);
%       
%       randomly pick k distinct numbers from a distribution p over 1:n.
%       Here, p should be a vector of length n, and it has sum(p) = 1.
%       

%   History
%   -------
%       - Created by Dahua Lin, on SEp 27, 2010
%

%% verify input

if isscalar(n)    
    use_p = false;
    
    k = double(k);
    n = double(n);
    
elseif isvector(n)
    
    p = n;
    n = numel(p);
    if size(p, 1) > 1
        p = p.';
    end    
    use_p = true;
    
    k = double(k);
    
else
    error('randpick:invalidarg', ...
        'The first argument should be either a scalar n or a vector p.');
end

if nargin < 3
    rstream = [];
end


%% main

if k > 0 && k < n
    
    if ~use_p
        if k == 1
            r = rp1(rstream, n);      
        
        elseif k == 2
            r = rp2(rstream, n);
  
        else              
            if n > 4 * k
                r = rpick_inc(rstream, n, k);
            elseif n > 4 * (n-k)
                r0 = rpick_inc(rstream, n, k);
                r = 1 : n;
                r(r0) = [];
            else
                r = rpick_sf(rstream, n, k);
            end

        end
        
    else % use_p
        
        if k == 1
            r = rdraw_p(rstream, p, 1);                        
        else
            r = rpick_inc_p(rstream, p, k);            
        end
        
    end
                    
elseif k == 0    
    r = [];
    
elseif k == n    
    r = 1 : n;
        
else
    error('randpick:invalidarg', ...
        'k is invalid, which should be with in the range [0, n].');    
end

    

%% core functions


function r = rpick_inc(rstream, n, k)
% use incremental sampling

r = uniqvalues(rdraw(rstream, n, k));
nr = numel(r);

while nr < k                
    s_new = uniqvalues(rdraw(rstream, n-nr, k-nr));
    r = rpick_merge_new(r, s_new);
    nr = numel(r);
end


function r = rpick_sf(rstream, n, k)
% use random-shuffling

if isempty(rstream)
    x = rand(1, n);
else
    x = rand(rstream, 1, n);
end

[sx, si] = sort(x); %#ok<ASGLU>
r = sort(si(1:k));


function r = rpick_inc_p(rstream, p, k)
% use incremental sampling with non-uniform p

r = uniqvalues(rdraw_p(rstream, p, k));
nr = numel(r);

while nr < k
    cp = p;
    cp(r) = [];
    cp = cp * (1 / sum(cp));
    
    s_new = uniqvalues(rdraw_p(rstream, cp, k-nr));
    r = rpick_merge_new(r, s_new);
    nr = numel(r);
end


%% auxiliary subfunction

function r = rp1(rstream, n)

if isempty(rstream)
    r = randi(n);
else
    r = randi(rstream, n);
end

function r = rp2(rstream, n)

if isempty(rstream)
    r = randi(n, 1, 2);
else
    r = randi(rstream, n, 1, 2);
end

if r(1) == r(2)  % re-draw r(2) from [1,2,...,r-1,n,r+1,...,n-1]
    r(2) = rp1(rstream, n-1);
    if r(2) == r(1)
        r(2) = n;
    end
end

if r(1) > r(2)
    r = [r(2) r(1)];
end


function r = rdraw(rstream, n, k)

if k == 1
    if isempty(rstream)
        r = randi(n);
    else
        r = randi(rstream, n);
    end
else
    if isempty(rstream)
        r = randi(n, 1, k);
    else
        r = randi(rstream, n, 1, k);
    end
end


function r = rdraw_p(rstream, p, k)
% F: cumulative distr. function of p

if isempty(rstream)
    x = rand(1, k);
else
    x = rand(rstream, 1, k);
end

F = cumsum(p);
r = pieces(x, F) + 1;


