function r = calc_sumprod(dim, n, a, b)
% sum the corresponding products along the direction dim
%
%   expand a or b to n rows/columns along the specified dimension
%

na = size(a, dim);
nb = size(b, dim);

if na == 1
    if nb == 1
        r = a .* b;
        if n > 1
            r = r * n;
        end
    else 
        r = a .* sum(b, dim);
    end
else
    if nb == 1
        r = sum(a, dim) .* b;
    else
        if dim == 1
            if size(a, 2) == 1
                r = a' * b;
            elseif size(b, 2) == 1
                r = b' * a;
            else
                r = sum(a .* b, 1);
            end
        else
            if size(a, 1) == 1
                r = b * a';
            elseif size(b, 1) == 1
                r = a * b';
            else
                r = sum(a .* b, 2);
            end            
        end
    end
end

        


