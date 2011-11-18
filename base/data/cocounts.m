function C = cocounts(siz, varargin)
% Counts the numbers of different integer combinations
%
%   C = cocounts(m, vs);
%       returns a vector of size m x 1, where C(k) being the number of 
%       occurrences of k in vs.
%
%   C = cocounts([m, n], vs1, vs2);
%       returns a matrix of size m x n, where 
%       C(i, j) = sum(vs1 == i & vs2 == j), which counts the occurrences
%       of the case where the value in vs1 is i and the corresponding
%       value in vs2 is j.
%
%   C = cocounts([n1, n2, n3, ...], vs1, vs2, vs3, ...);
%       returns an array C of size n1 x n2 x n3 x ... in a similar manner.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 10, 2010
%

%% verify input arguments

if ~(isnumeric(siz) && ndims(siz) == 2 && size(siz,1) == 1)
    error('cocounts:invalidarg', 'siz should be a row vector.');
end

d = numel(siz);

if numel(varargin) ~= d
    error('cocounts:invalidarg', 'The number of value vectors is invalid.');
end


%% main

if d == 1
    vs = varargin{1};    
else
    vs = sub2ind(siz, varargin{:});
end

N = prod(siz);
C = intcount(N, vs);

if d == 1
    C = C(:);
else
    C = reshape(C, siz);
end

