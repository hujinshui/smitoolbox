% Count the number of occurrence of integers
%
%   c = intcount([v0, v1], V);
%       counts the number of integers in [v0, v1] that occurs in V.
%       
%       In the input, V is an array containing integer values (whose
%       type can be double or single nonetheless). 
%       In the output, c is a row vector of length v1 - v0 + 1.
%       In particular, c(i) counts the number of times that the value
%       v0 + (i-1) occurs in V.
%
%   c = intcount(n, V);
%       This is equivalent to intcount([1, n], V);
%

%   History
%   -------
%       - Created by Dahua Lin, on May 26, 2010
%       - Modified by Dahua Lin, on Jun 6, 2010
%           - Use pure C++ mex implementation
%
