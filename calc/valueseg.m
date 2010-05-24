function [sp, ep, v] = valueseg(x)
%VALUESEG divides a vector of numeric values into constant-value segments
%
%   [sp, ep] = valueseg(x);
%       divides the numeric vector x into ranges, such that the elements
%       in each range have the same value.
%
%       In the output, sp and ep are vectors of equal sizes, whose length
%       equals the number of segments. They respectively give the start
%       and end index of the segments.
%
%   [sp, ep, v] = valueseg(x);
%       this additionally returns the value shared by the elements of 
%       each segment. 
%

%   History
%       - Created by Dahua Lin, on May 30, 2008
%

%% main

if isempty(x)
    sp = [];
    ep = [];
    return
end

assert(isnumeric(x) && isvector(x), 'valueseg:invalidarg', ...
    'x should be a numeric vector.');

n = numel(x);

if n == 1
    sp = 1;
    ep = 1;
    
else
    dp = find(diff(x));
    
    if size(x, 1) == 1      % row vector
        sp = [1, 1 + dp];
        ep = [dp, n];
    else                    % column vector
        sp = [1; 1 + dp];
        ep = [dp; n];
    end
end

if nargout >= 3
    v = x(sp);
end


