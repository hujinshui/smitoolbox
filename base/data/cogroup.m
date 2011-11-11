function G = cogroup(siz, varargin)
% Groups the indices for different integer combinations
%
%   G = cogroup(m, vs);
%       returns a cell vector of size m x 1, where G{k} is a vector of
%       indices, such that all(vs(G{k}) == k).
%
%   C = cogroup([m, n], vs1, vs2);
%       returns a matrix of size m x n, where all(vs1(G{i,j}) == i) and
%       all(vs2(G{i,j}) == j).
%
%   C = cogroup([n1, n2, n3, ...], vs1, vs2, vs3, ...);
%       returns an array C of size n1 x n2 x n3 x ... in a similar manner.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 10, 2010
%

%% verify input arguments

if ~(isnumeric(siz) && ndims(siz) == 2 && size(siz,1) == 1)
    error('cogroup:invalidarg', 'siz should be a row vector.');
end

d = numel(siz);

if numel(varargin) ~= d
    error('cogroup:invalidarg', 'The number of value vectors is invalid.');
end


%% main

if d == 1
    vs = varargin{1};    
else
    vs = sub2ind(siz, varargin{:});
end

N = prod(siz);
G = intgroup(N, vs);

if d == 1
    G = G(:);
else
    G = reshape(G, siz);
end

