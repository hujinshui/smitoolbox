function [sp, ep, v] = valueseg(x)
%VALUESEG divides a vector of numeric values into constant-value segments
%
%   [sp, ep] = valueseg(x);
%       divides the numeric/logical/char vector x into ranges, such that 
%       the elements in each range have the same value.
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
%       - Modified by Dahua Lin, on May 26, 2010
%           - based on C++-mex implementation
%

%% verify input

if ~isvector(x) || isempty(x) || issparse(x)
    error('valueseg:invalidarg', 'x should be a non-empty full vector.');
end

%% main

[sp, ep] = valueseg_cimp(x);

if nargout >= 3
    v = x(sp);
end


