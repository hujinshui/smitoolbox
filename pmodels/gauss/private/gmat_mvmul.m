function R = gmat_mvmul(cf, C, V)
% compute matrix-vector multiplication
%

switch cf
    case 's'
        n = size(C, 2);
        if n == 1
            R = C * V;
        else
            R = bsxfun(@times, C, V);
        end
        
    case 'd'
        n = size(C, 2);
        if n == 1
            R = bsxfun(@times, C, V);
        else
            R = C .* V;
        end
        
    case 'f'
        n = size(C, 3);
        if n == 1
            R = C * V;
        else
            R = zeros(size(C,1), n, fresult_class(C, V));
            for i = 1 : n
                R(:,i) = C(:,:,i) * V(:,i);
            end
        end
end

