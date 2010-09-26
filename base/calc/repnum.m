% Generates a vector by repeating numbers
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

%   History
%   -------
%       - Created by Dahua Lin, on Nov 3, 2009
%       - Modified by Dahua Lin, on Jun 6, 2010
%           - use C++ mex to improve efficiency
%

