function R = gmat_choltrans(cf, C, X)
% Compute chol(C) * X
%
% There should be only one covariance packed in C
%

switch cf
    case 's'
        R = sqrt(C) * X;
        
    case 'd'
        d = size(C, 1);
        if d == 1
            R = sqrt(C) * X;
        else
            R = bsxfun(@times, sqrt(C), X);
        end
        
    case 'f'
        d = size(C, 1);
        if d == 1
            R = sqrt(C) * X;
        else
            R = chol(C) * X;
        end
            
end