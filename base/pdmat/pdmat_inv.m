function R = pdmat_inv(S)
% Compute the inverse of the positive definite matrices
%
%   R = pdmat_inv(S);
%       return the inverse of the positive definite matrices in S.
%       The inverse is represented in the same form as in S.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% main

R.tag = S.tag;
R.ty = S.ty;
R.d = S.d;
R.n = S.n;

switch S.ty
    case {'s', 'd'}
        R.v = 1 ./ S.v;
        
    case 'f'
        if S.d == 1
            R.v = 1 ./ S.v;
            
        elseif S.d == 2
            R.v = inv2x2(S.v);
        
        else
            if S.n == 1
                R.v = inv(S.v);
            else
                v = S.v;
                r = zeros(size(S.v), class(S.v));
                for i = 1 : S.n
                    r(:,:,i) = inv(v(:,:,i));
                end
                R.v = r;
            end
        end        
end

