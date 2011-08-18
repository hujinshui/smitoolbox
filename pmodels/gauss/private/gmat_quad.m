function R = gmat_quad(cf, C, X, Y)
% compute quadratic form w.r.t. Gaussian matrix
%
%   R would be of size m x n, where
%   m - # matrices packed in C
%   n - # columns in X or Y
%

d = size(X, 1);

switch cf
    case 's'
        if d == 1
            v = X .* Y;
        else
            v = dot(X, Y, 1);
        end
        if isscalar(C)
            R = C * v;
        else
            R = bsxfun(@times, C.', v);
        end            
        
    case 'd'
        R = C' * (X .* Y);
        
    case 'f'
        m = size(C, 3);
        if m == 1
            if d == 1
                R = X .* (C * Y);
            else
                R = dot(X, C * Y, 1);
            end                
        else
            R = zeros(m, size(X,2), fresult_class(C, X, Y));
            if d == 1
                for i = 1 : m
                    R(i, :) = X .* (C(:,:,i) .* Y);
                end
            else
                for i = 1 : m
                    R(i, :) = dot(X, C(:,:,i) * Y, 1);
                end
            end
        end
end
