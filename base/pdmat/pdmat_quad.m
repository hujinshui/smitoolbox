function R = pdmat_quad(S, X, Y)
% computes quadratic form with respect to positive definite matrices
%
%   R = pdmat_quad(S, X);
%   R = pdmat_quad(S, X, Y);
%       
%       Supposes there are m matrices contained in S, and X and Y
%       are matrices of the same size S.d x n.
%
%       Then R is a matrix of size m x n, with R(i, j) being equal to
%       x(:,i)' * Ai * y(:, i), where Ai is the i-th matrix in S.
%
%       When Y is omitted, it is assumed to be equal to X.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

if nargin < 3
    Y = X;
end

if ~(isfloat(X) && ndims(X) == 2)
    error('pdmat_quad:invalidarg', 'X should be a numeric matrix.');
end

if ~(isfloat(Y) && ndims(Y) == 2)
    error('pdmat_quad:invalidarg', 'Y should be a numeric matrix.');
end

d = S.d;
if ~(size(X,1) == d && size(Y,1) == d && size(X,2) == size(Y,2))
    error('pdmat_quad:invalidarg', ...
        'X and Y should be equal-size matrices with S,d rows.');
end


%% main

m = S.n;
n = size(X, 2);

switch S.ty
    case 's'
        if d == 1
            xyd = X .* Y;
        else
            xyd = dot(X, Y, 1);
        end
        R = S.v' * xyd;
        
    case 'd'
        R = S.v' * (X .* Y);
        
    case 'f'
        v = S.v;
        if d == 1
            if m > 1
                v = reshape(v, m, 1);
            end
            R = v * (X .* Y);
        else
            if m == 1
                R = dot(X, v * Y, 1);
            else
                r1 = dot(X, v(:,:,1) * Y, 1);
                R = zeros(m, n, class(r1));
                R(1,:) = r1;
                for k = 2 : m
                    R(k,:) = dot(X, v(:,:,k) * Y, 1);
                end
            end
        end
end


