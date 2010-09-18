%VALUESEG divides a vector of numeric values into constant-value segments
%
%   sp = valueseg(x);
%       divides the numeric/logical/char vector x into segments, such that 
%       the elements in each segments have the same value.
%
%       In the output, sp is a row vector of size m, where m is the number
%       of segments. In particular, sp(i) is the starting index of the
%       i-th value segment.
%

%   History
%       - Created by Dahua Lin, on May 30, 2008
%       - Modified by Dahua Lin, on May 26, 2010
%           - based on C++-mex implementation
%       - Modified by Dahua Lin, on June 6, 2010
%           - change to pure C++-mex implementation
%



