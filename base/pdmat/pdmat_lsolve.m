function Y = pdmat_lsolve(S, X)
% Solves linear equation solution w.r.t positive definite matrices
%
%   Y = pdmat_lsolve(S, X);
%       
%       Returns the result Y as pdmat_mvmul(inv(S), X), with the
%       calculation implemented in a more efficient way.
%

% Created by Dahua Lin, on Aug 25, 2010
%


%% verify input

d = S.d;
n = S.n;

if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == d)
    error('pdmat_lsolve:invalidarg', ...
        'X should be a numeric matrix with size(X,1) == d.');
end

%% main

if n == 1
    
    switch S.ty
        case 's'
            Y = X / S.v;
        case 'd'
            if d == 1
                Y = X / S.v;
                
            elseif size(X, 2) == 1
                Y = X ./ S.v;
                
            else
                Y = bsxfun(@times, 1 ./ S.v, X);
            end
        case 'f'
            Y = S.v \ X;
    end
    
else % n > 1
    
    if size(X, 2) ~= n
        error('pdmat_mvmul:invalidarg', ...
            '#columns in X is not consistent with S.n.');
    end
    
    switch S.ty
        case 's'
            if d == 1
                Y = X ./ S.v;
            else
                Y = bsxfun(@times, 1 ./ S.v, X);
            end
            
        case 'd'
            Y = X ./ S.v;
            
        case 'f'
            if d == 1
                Y = X ./ reshape(S.v, 1, n);
            else
                v = S.v;
                y1 = v(:,:,1) \ X(:,1);
                Y = zeros(d, n, class(y1));
                Y(:, 1) = y1;
                for i = 2 : n
                    Y(:, i) = v(:,:,i) \ X(:,i);
                end
            end
    end
end