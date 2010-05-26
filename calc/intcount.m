function c = intcount(rgn, V)
% Count the number of occurrence of integers
%
%   c = intcount([v0, v1], v);
%       counts the number of integers in [v0, v1] that occurs in V.
%       
%       In the input, V is an array containing integer values (whose
%       type can be double or single nonetheless). 
%       In the output, c is a row vector of length v1 - v0 + 1.
%       In particular, c(i) counts the number of times that the value
%       v0 + (i-1) occurs in V.
%

% Created by Dahua Lin, on May 26, 2010
%

%% verify input arguments

if ~(isnumeric(rgn) && numel(rgn) == 2)
    error('intcount:invalidarg', ...
        'The value range should be given in form of a pair of numbers.');
end

if ~(isnumeric(V))
    error('intcount:invalidarg', 'V should be a numeric array.');
end

%% main

if ~isa(V, 'int32')
    V = int32(V);
end

c = intcount_cimp(double(rgn), V);

