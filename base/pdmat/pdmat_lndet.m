function r = pdmat_lndet(S)
% Compute log-determinant of positive definite matrices
%
%   r = pdmat_lndet(S);
%       computes the log-determinant of positive definite matrices
%       contained in S.
%
%       r is a scalar (if S.n == 1) or a row vector of size 1 x S.n,
%       with r(i) corresponding to the i-th matrix in S.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% main

d = S.d;

switch S.ty
    case 's'
        r = d * log(S.v);
        
    case 'd'
        if d == 1
            r = log(S.v);
        else            
            r = sum(log(S.v), 1);
        end
                
    case 'f'
        if d == 1
            r = log(reshape(S.v, 1, S.n));
            
        elseif d == 2
            r = log(det2x2(S.v));
        
        else
            n = S.n;
            v = S.v;
            if n == 1
                r = lndet(v);
            else
                r = zeros(1, n, class(v));
                for i = 1 : n
                    r(i) = lndet(v(:,:,i));
                end
            end
        end        
end

