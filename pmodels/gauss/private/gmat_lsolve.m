function R = gmat_lsolve(cf, C, V)
% Solve linear equations
%

switch cf
    case 's'
        n = size(C, 2);
        if n == 1
            R = (1/C) * V;
        else
            R = bsxfun(@times, 1 ./ C, V);
        end
        
    case 'd'
        n = size(C, 2);
        if n == 1
            R = bsxfun(@times, 1 ./ C, V);
        else
            R = V ./ C;
        end
        
    case 'f'
        n = size(C, 3);
        if n == 1
            R = C \ V;
        else
            R = zeros(size(C,1), n, fresult_class(C, V));
            for i = 1 : n
                R(:,i) = C(:,:,i) \ V(:,i);
            end
        end
end