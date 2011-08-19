function R = gmat_fullform(cf, d, C)
% Get the full covariance matrix of a given representation
%
%   R = gmat_fullform(cf, d, C);
%

switch cf
    case 's'
        n = numel(C);
        if n == 1
            R = diag(C * ones(1, d));
        else            
            R = zeros(d * d, n, class(C));
            R(1:(d+1):d*d, :) = C(ones(d, 1), :);
            R = reshape(R, [d, d, n]);
        end
        
    case 'd'
        n = size(C, 2);
        if n == 1
            R = diag(C);
        else
            R = zeros(d * d, n, class(C));
            R(1:(d+1):d*d, :) = C;
            R = reshape(R, [d, d, n]);
        end
        
    case 'f'
        R = C;
        
end