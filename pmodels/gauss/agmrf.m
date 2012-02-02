function G = agmrf(W, a, y)
% Constructs an attractive Gaussian MRF
%
%   An attractive Gaussian MRF is defined to be a Gaussian model, as
%
%       p(x) = exp(E(x)) / Z,
%
%   with
%
%       E(x) = 1/2 * sum_{(i,j) \in E} W(i,j) (x(i) - x(j))^2
%            + 1/2 * a (x(i) - y(i))^2
%
%   This can be written in a matrix form as
%
%       E(x) = 1/2 * (x' * L * x + (x-y)' * Da * (x-y) ).
%
%   G = agmrf(W, a);
%   G = agmrf(W, a, y);
%
%       Constructs an attractive Gaussian MRF as formulated above. 
%       
%       Input arguments:
%       - W:        The weighted adjacency matrix [n x n]
%       - a:        The weights that link x to y, which can be either a
%                   scalar or a vector of length n.
%       - y:        The y values, which can be a vector of length n.
%                   If y is omitted, it is assumed to be zero.
%
%       In the output, G is a gaussd struct with G.ty == 'c' (using
%       canonical parameters).
%

% Created by Dahua Lin, on Dec 9, 2011
%

%% verify inputs

if ~(isfloat(W) && ndims(W) == 2 && isreal(W))
    error('agmrf:invalidarg', 'W should be a real matrix.');
end
n = size(W, 1);

if ~(isfloat(a) && isreal(a) && ...
        (isscalar(a) || (isvector(a) && length(a) == n)))
    error('agmrf:invalidarg', 'a should be a real vector of length n.');
end

if nargin < 3
    y = 0;
else
    if ~(isfloat(y) && isreal(y) && isvector(y) && length(y) == n)
        error('agmrf:invalidarg', 'y should be a real vector of length n.');
    end
end
            
%% main

J = laplacemat(W, a);

if isequal(y, 0) || isequal(a, 0)
    h = 0;
else
    if size(y, 2) > 1; y = y.'; end
    if size(a, 2) > 1; a = a.'; end
    
    if isequal(a, 1)
        h = y;
    else
        h = y .* a;
    end
end

G = gaussd('c', h, J);
        

