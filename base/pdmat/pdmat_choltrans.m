function Y = pdmat_choltrans(S, X)
% Transform vectors using the result of Cholesky factorization
%
%   Y = pdmat_choltrans(S, X);
%
%       Here, S should contain only one matrix (denoted by A). 
%       Let L be a lower-triangular matrix such that A = L * L'.
%       Then this function returns Y = L * X, which is computed
%       in an efficient manner, depending on the form of S.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

if S.n ~= 1
    error('pdmat_pwquad:invalidarg', 'S must be comprised of one matrix.');
end

d = S.d;
if ~(isfloat(X) && ndims(X) == 2 && size(X,1) == d)
    error('pdmat_pwquad:invalidarg', ...
        'X should be a numeric matrix with size(X,1) == d.');
end

%% main

switch S.ty
    case 's'
        Y = sqrt(S.v) * X;
        
    case 'd'
        if d == 1 
            Y = sqrt(S.v) * X;            
        elseif size(X,2) == 1
            Y = sqrt(S.v) .* X;
        else
            Y = bsxfun(@times, sqrt(S.v), X);
        end
        
    case 'f'
        if d == 1
            Y = sqrt(S.v) * X;
        elseif d == 2
            Y = chol2x2(S.v) * X;
        else
            Y = chol(S.v, 'lower') * X;
        end
end

