function Y = pdmat_mvmul(S, X)
% Compute product of positive-definite matrix and vectors
%
%   Y = pdmat_mvmul(S, X)
%       compute the product of the positive definite matrix (matrices) 
%       and the columns of X.
%
%       Several cases are supports (let ns = S.n and nx = size(X, 2))
%       - ns == 1 && nx == 1: 
%         times the matrix represented by S with the vector x.
%         size(Y) = [d, 1]
%         
%       - ns == 1 && nx > 1:
%         times the matrix represented by S with all columns of x.
%         size(Y) = [d, nx]
%
%       - ns == nx > 1:
%         times the i-th matrix contained in S with X(:,i).
%         size(Y) = [d, nx]
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

d = S.d;
n = S.n;

if ~(isfloat(X) && ndims(X) == 2 && size(X, 1) == d)
    error('pdmat_mvmul:invalidarg', ...
        'X should be a numeric matrix with size(X,1) == d.');
end

%% main

if n == 1
    
    switch S.ty
        case {'s', 'f'}
            Y = S.v * X;
        case 'd'
            if d == 1
                Y = S.v * X;
                
            elseif size(X, 2) == 1
                Y = S.v .* X;
                
            else
                Y = bsxfun(@times, S.v, X);
            end
    end
    
else % n > 1
    
    if size(X, 2) ~= n
        error('pdmat_mvmul:invalidarg', ...
            '#columns in X is not consistent with S.n.');
    end
    
    switch S.ty
        case 's'
            if d == 1
                Y = S.v .* X;
            else
                Y = bsxfun(@times, S.v, X);
            end
            
        case 'd'
            Y = S.v .* X;
            
        case 'f'
            if d == 1
                Y = reshape(S.v, 1, n) .* X;
            else
                v = S.v;
                y1 = v(:,:,1) * X(:,1);
                Y = zeros(d, n, class(y1));
                Y(:, 1) = y1;
                for i = 2 : n
                    Y(:, i) = v(:,:,i) * X(:,i);
                end
            end
    end
end


