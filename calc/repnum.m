function x = repnum(ns)
%Generates a vector by repeating numbers
%
%   x = repnum(ns);
%       generates a vector x, whose length is sum(ns), by repeating i for
%       ns(i) times.
%
%       For example, repnum([3 3 4]) generates a vector of length 10 as
%       [1 1 1 2 2 2 3 3 3 3].
%
%       The result of repnum can be used as indices for generating 
%       other arrays with repeated values. 
%
%       Example
%       -------
%           x = [0.1 0.2 0.4]
%           x(repnum([2 2 3])) is [0.1 0.1 0.2 0.2 0.4 0.4 0.4].
%

% Created by Dahua Lin, on Nov 3, 2009
%

%% parse and verify input arguments

assert(isnumeric(ns) && isvector(ns), 'repnum:invalidarg', ...
    'ns should be a numeric vector.');

%% main

if size(ns, 1) > 1
    tt = true;
    ns = ns.';
else
    tt = false;
end

if ~all(ns)
    v = find(ns);
    ns = ns(v);
    uv = true;
else
    uv = false;
end

N = sum(ns);
x = zeros(1, N);

x([1, 1 + cumsum(ns(1:end-1))]) = 1;

x = cumsum(x);
if uv
    x = v(x);
end

if tt
    x = x.';
end

