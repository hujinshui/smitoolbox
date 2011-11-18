function V = pdmat_diag(S)
% Get the diagonal entries of the matrices
%
%   V = pdmat_diag(S);
%
%       returns the diagonal entries of each matrices packed in S
%       as columns of V. 
%

% Created by Dahua Lin, on Sep 1, 2010
%

switch S.ty
    case 's'
        if S.d == 1
            V = S.v;
        else
            V = S.v(ones(S.d, 1), :);
        end
    case 'd'
        V = S.v;
    case 'f'
        if S.d == 1
            if S.n == 1
                V = S.v;
            else
                V = reshape(S.v, [1, S.n]);
            end
        else
            if S.n == 1
                V = diag(S.v);
            else
                d = S.d;
                di = bsxfun(@plus, (1:(d+1):(d*d)).', (0:S.n-1) * (d*d));
                V = S.v(di);
            end
        end
end


