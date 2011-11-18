function R = pdmat_pwquad(S, X, Y)
% computes pairwise quadratic forms w.r.t. positive definite matrices
%
%   R = pdmat_pwquad(S, X);
%   R = pdmat_pwquad(S, X, Y);
%       
%       Here, S must be comprised of exactly one matrix (S.n == 1).
%       
%       Suppose the size of X is S.d x m and that of Y is S.d x n. Then
%       R will be an m x n matrix, with R(i, j) = X(:,i)' * A * Y(:,i).
%       Here, A denotes the matrix represented by S.
%
%       When Y is omitted, it is assumed to be equal to X. In this case,
%       R is guaranteed to be exactly symmetric.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

if nargin < 3
    is_sym = 1;
else
    is_sym = 0;
end

if S.n ~= 1
    error('pdmat_pwquad:invalidarg', 'S must be comprised of one matrix.');
end

d = S.d;
if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == d)
    error('pdmat_pwquad:invalidarg', ...
        'X should be a numeric matrix with size(X,1) == d.');
end

if ~is_sym
    if ~(isfloat(Y) && ndims(Y) == 2 && size(Y,1) == d)
        error('pdmat_pwquad:invalidarg', ...
            'Y should be a numeric matrix with size(Y,1) == d.');
    end
end

%% main

switch S.ty
    case 's'
        if is_sym
            R = S.v * (X' * X);
        else
            R = S.v * (X' * Y);
        end            
        
    case 'd'
        if is_sym
            R = X' * bsxfun(@times, S.v, X);
            R = (R + R') * 0.5;
        else
            R = X' * bsxfun(@times, S.v, Y);
        end

    case 'f'
        if is_sym
            R = X' * (S.v * X);
            R = (R + R') * 0.5;
        else
            R = X' * (S.v * Y);
        end
end


