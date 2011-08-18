function R = gmat_choltrans(cf, C, X)
% Compute chol(C) * X
%
% There should be only one covariance packed in C
%

switch cf
    case 's'
        if C == 1
            R = X;
        else
            R = sqrt(C) * X;
        end
        
    case 'd'
        d = size(C, 1);
        if d == 1
            if C == 1
                R = X;
            else
                R = sqrt(C) * X;
            end
        else
            R = bsxfun(@times, sqrt(C), X);
        end
        
    case 'f'
        d = size(C, 1);
        if d == 1
            R = sqrt(C) * X;
        elseif d == 2
            R = chol2x2(C) * X;
        else
            R = chol(C, 'lower') * X;
        end
            
end