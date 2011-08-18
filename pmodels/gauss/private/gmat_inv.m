function R = gmat_inv(cf, C)
% compute inverse of Gaussian matrix
%

switch cf
    case 's'
        R = 1 ./ C;
    case 'd'
        R = 1 ./ C;
    case 'f'
        n = size(C, 3);
        d = size(C, 1);
        
        if d == 1
            R = 1 ./ C;
        elseif d == 2
            R = inv2x2(C);
        else
            if n == 1
                R = inv(C);
            else
                R = zeros(size(C), class(C));
                for i = 1 : n
                    R(:,:,i) = inv(C(:,:,i));
                end
            end
        end
end
