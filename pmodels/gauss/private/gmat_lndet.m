function ldc = gmat_lndet(cf, d, C)
% Compute the log of determinant
%

switch cf
    case 's'
        ldc = d * log(C);
    case 'd'
        ldc = sum(log(C), 1);
    case 'f'
        n = size(C, 3);
        if n == 1
            ldc = lndet(C);
        else
            ldc = zeros(1, n, class(C));
            for i = 1 : n
                ldc(i) = lndet(C(:,:,i));                
            end
        end        
end