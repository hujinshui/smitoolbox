function dv = pdmat_dot(A, B)
% Computes the inner product between positive definite matrices
%
%   dv = pdmat_dot(A);
%   dv = pdmat_dot(A, B);
%
%       Computes the inner product (tr(A' * B)) between A and B.
%
%       The following cases are supported:
%       - A.n == 1 && B.n == 1:     return dv as a scalar
%       - A.n == 1 && B.n = n > 1:  return dv as a 1 x n row vector
%       - A.n = n > 1 && B.n == 1:  return dv as a 1 x n row vector
%       - A.n == B.n = n > 1:       return dv as a 1 x n row vector.
%
%       If B is omitted, it sets B = A.
%

% Created by Dahua Lin, on Sep 2, 2011
%

%% verify input

d = A.d;

if nargin < 2
    B = A;
else
    if B.d ~= d
        error('pdmat_dot:invalidarg', ...
            'The dimensions of A and B are inconsistent.');
    end
    
    if ~(A.n == B.n || A.n == 1 || B.n == 1)
        error('pdmat_dot:invalidarg', ...
            'The numbers of matrices in A and B are inconsistent.');
    end            
end

%% main

switch A.ty
    case 's'
        switch B.ty
            case 's'
                dv = dot_ss(d, A, B);
            case 'd'
                dv = dot_sd(d, A, B);
            case 'f'
                dv = dot_sf(d, A, B);
        end
        
    case 'd'
        switch B.ty
            case 's'
                dv = dot_sd(d, B, A);
            case 'd'
                dv = dot_dd(d, A, B);
            case 'f'
                dv = dot_df(d, A, B);
        end
        
    case 'f'
        switch B.ty
            case 's'
                dv = dot_sf(d, B, A);
            case 'd'
                dv = dot_df(d, B, A);
            case 'f'
                dv = dot_ff(d, A, B);
        end
        
end


%% Core computation routines

function dv = dot_ss(d, A, B)

dv = A.v .* B.v;
if d > 1
    dv = dv * d;
end

function dv = dot_sd(d, A, B)

if d == 1
    dv = A.v .* B.v;
else
    dv = A.v .* sum(B.v, 1);
end

function dv = dot_dd(d, A, B)

if d == 1
    dv = A.v .* B.v;
else
    dv = col_dots(A.v, B.v);
end

function dv = dot_sf(d, A, B)

if d == 1
    if B.n == 1
        bv = B.v;
    else
        bv = reshape(B.v, 1, B.n);
    end
    dv = A.v .* bv;
else
    bdiagv = get_diagv(B.v);
    dv = A.v .* sum(bdiagv, 1);
end
        
function dv = dot_df(d, A, B)

if d == 1
    if B.n == 1
        bv = B.v;
    else
        bv = reshape(B.v, 1, B.n);
    end
    dv = A.v .* bv;
else
    bdiagv = get_diagv(B.v);
    dv = col_dots(A.v, bdiagv);
end

function dv = dot_ff(d, A, B)

if A.n == 1
    av = A.v(:);
else
    av = reshape(A.v, d * d, A.n);
end

if B.n == 1
    bv = B.v(:);
else
    bv = reshape(B.v, d * d, B.n);
end

if d == 1
    dv = av .* bv;
else
    dv = col_dots(av, bv);
end

        
%% Auxiliary functions

function V = get_diagv(A)
% get diagonal entries of fullform matrix or matrices

n = size(A, 3);
if n == 1
    V = diag(A);
else
    d = size(A, 1);
    I = bsxfun(@plus, (1:(d+1):(d*d))', (0:n-1) * (d*d));
    V = A(I);
end

function v = col_dots(X, Y)

if size(X, 2) == 1
    v = X' * Y;
elseif size(Y, 2) == 1
    v = Y' * X;
else
    v = dot(X, Y, 1);
end






