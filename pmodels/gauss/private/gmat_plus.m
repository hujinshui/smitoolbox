function [Cr, cfr] = gmat_plus(C1, cf1, C2, cf2)
% Compute the sum of two covariance
%
%   [Cr, cfr] = gmat_plus(C1, cf1, C2, cf2);
%
%       C1:     the first input covariance of form cf1
%       C2:     the second input covariance of form cf2
%
%       C2:     the resultant covariance of form cfr
%

%% main delegation

switch cf1
    case 's'
        switch cf2
            case 's'
                [Cr, cfr] = plus_ss(C1, C2);
            case 'd'
                [Cr, cfr] = plus_ds(C2, C1);
            case 'f'
                [Cr, cfr] = plus_fs(C2, C1);
        end
    case 'd'
        switch cf2
            case 's'
                [Cr, cfr] = plus_ds(C1, C2);
            case 'd'
                [Cr, cfr] = plus_dd(C1, C2);
            case 'f'
                [Cr, cfr] = plus_fd(C2, C1);
        end
    case 'f'
        switch cf2
            case 's'
                [Cr, cfr] = plus_fs(C1, C2);
            case 'd'
                [Cr, cfr] = plus_fd(C1, C2);
            case 'f'
                [Cr, cfr] = plus_ff(C1, C2);
        end
end

%% core routines


function [Cr, cfr] = plus_ss(C1, C2)

cfr = 's';
Cr = C1 + C2;

function [Cr, cfr] = plus_dd(C1, C2)

cfr = 'd';
if size(C1, 2) == size(C2, 2)
    Cr = C1 + C2;
else
    Cr = bsxfun(@plus, C1, C2);
end

function [Cr, cfr] = plus_ff(C1, C2)

cfr = 'f';
if size(C1, 3) == size(C2, 3)
    Cr = C1 + C2;
else
    Cr = bsxfun(@plus, C1, C2);
end

function [Cr, cfr] = plus_ds(C1, C2)

cfr = 'd';
d = size(C1, 1);
if d == 1
    Cr = C1 + C2;
else
    if isscalar(C2)
        Cr = C1 + C2;
    else
        Cr = bsxfun(@plus, C1, C2);
    end
end

function [Cr, cfr] = plus_fs(C1, C2)

cfr = 'f';
d = size(C1, 1);
n1 = size(C1, 3);
n2 = numel(C2);

if d == 1
    if n2 == 1
        Cr = C1 + C2;
    elseif n1 == 1
        Cr = reshape(C1 + C2, [1, 1, n2]);
    else
        Cr = C1 + reshape(C2, [1, 1, n2]);
    end
else
    if n1 == 1 && n2 == 1
        Cr = adddiag(C1, C2);
    else
        n = max(n1, n2);
        Cr = reshape(C1, d*d, n1);
        if n1 == 1
            Cr = repmat(Cr, 1, n);
        end
        if n2 == 1
            Cr(1:(d+1):d*d, :) = Cr(1:(d+1):d*d, :) + C2;
        else
            Cr(1:(d+1):d*d, :) = bsxfun(@plus, Cr(1:(d+1):d*d, :), C2);
        end
        
        Cr = reshape(Cr, [d, d, n]);  
    end
end
        

function [Cr, cfr] = plus_fd(C1, C2)

cfr = 'f';
d = size(C1, 1);
n1 = size(C1, 3);
n2 = size(C2, 2);

if d == 1
    if n2 == 1
        Cr = C1 + C2;
    elseif n1 == 1
        Cr = reshape(C1 + C2, [1, 1, n2]);
    else
        Cr = C1 + reshape(C2, [1, 1, n2]);
    end
else
    if n1 == 1 && n2 == 1
        Cr = adddiag(C1, C2);
    else
        n = max(n1, n2);
        Cr = reshape(C1, d*d, n1);
        if n1 == 1
            Cr = repmat(Cr, 1, n);
        end
        if n2 == n
            Cr(1:(d+1):d*d, :) = Cr(1:(d+1):d*d, :) + C2;
        else
            Cr(1:(d+1):d*d, :) = bsxfun(@plus, Cr(1:(d+1):d*d, :), C2);
        end
        
        Cr = reshape(Cr, [d, d, n]);
    end
end




