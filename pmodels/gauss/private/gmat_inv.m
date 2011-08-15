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
        if n == 1
            R = inv(C);
        else
            R = zeros(size(C), class(C));
            for i = 1 : n
                C(:,:,i) = inv(R(:,:,i));
            end
        end
end
