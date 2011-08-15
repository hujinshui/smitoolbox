function R = gmat_quad(cf, C, X, Y)
% compute quadratic form w.r.t. Gaussian matrix
%
%   R would be of size m x n, where
%   m - # matrices packed in C
%   n - # columns in X or Y
%

switch cf
    case 's'
        v = dot(X, Y, 1);
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
            R = dot(X, C * Y, 1);
        else
            R = zeros(m, size(X,2), fresult_class(C, X, Y));
            for i = 1 : m
                R(:,:,i) = dot(X, C(:,:,i) * Y, 1);
            end
        end
end
