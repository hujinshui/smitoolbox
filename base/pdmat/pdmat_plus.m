function Sr = pdmat_plus(S1, S2, c1, c2)
% Compute the sum or linear combination of two matrices
%
%   Sr = pdmat_plus(S1, S2);
%       computes the sum of the matrix in S1 and that in S2, and returns
%       the resultant matrix in a pdmat struct.
%
%   Sr = pdmat_plus(S1, S2, c1, c2);
%       computes the linear combination of the matrix in S1 (scaled with
%       c1) and that in S2 (scaled with c2).
%
%   Several cases are supported, as below.
%   - S1.n == S2.n = 1 or any n > 0
%   - S1.n == 1 and S2.n > 0
%   - S1.n > 0 and S2.n == 1
%
%   S1 and S2 may or may not be in the same form, the results are in a
%   promoted form as max(S1.ty, S2.ty) under the order 's' < 'd' < 'f'.
%       '

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

n1 = S1.n;
n2 = S2.n;

d = S1.d;
if d ~= S2.d
    error('pdmat_plus:invalidarg', 'S1.d and S2.d are inconsistent.');
end

if ~(n1 == n2 || n1 == 1 || n2 == 1)
    error('pdmat_plus:invalidarg', ...
        'The number of matrices in S1 and that in S2 are incompatible.');
end

if nargin == 2
    c1 = 1;
    c2 = 1;
elseif nargin == 4
    if ~(isfloat(c1) && isscalar(c1) && isreal(c1))
        error('pdmat_plus:invalidarg', 'c1 must be a real scalar.');
    end
    if ~(isfloat(c2) && isscalar(c2) && isreal(c2))
        error('pdmat_plus:invalidarg', 'c2 must be a real scalar.');
    end
else
    error('The #arguments to pdmat_plus must be either 2 or 4.');
end


%% main

% extract & scale data

if c1 == 1
    v1 = S1.v;
else
    v1 = c1 * S1.v;
end

if c2 == 1
    v2 = S2.v;
else
    v2 = c2 * S2.v;
end

% do addition

ty = max_ty(S1.ty, S2.ty);

if n1 == n2 
    n = n1;
    vr = calc_sn(d, n, S1.ty, S2.ty, v1, v2);
           
else    
    if n1 < n2
        n = n2;
        vr = calc_dn(d, n, S1.ty, S2.ty, v1, v2);
    else
        n = n1;
        vr = calc_dn(d, n, S2.ty, S1.ty, v2, v1);
    end    
end

% make result

Sr.tag = S1.tag;
Sr.ty = ty;
Sr.d = d;
Sr.n = n;
Sr.v = vr;

    

%% sub-routines

function ty = max_ty(ty1, ty2)

if ty1 == 's' || (ty1 == 'd' && ty2 == 'f')
    ty = ty2;
else
    ty = ty1;
end     


function vr = calc_sn(d, n, ty1, ty2, v1, v2)
% pre-conditions:
%   S1.n == S2.n

if ty1 == ty2
    vr = v1 + v2;
    
else
    if ty1 == 's'
        if ty2 == 'd'
            vr = add_sd(d, n, v1, v2);
        elseif ty2 == 'f'
            vr = add_sf(d, n, v1, v2);
        end
        
    elseif ty1 == 'd'
        if ty2 == 's'
            vr = add_sd(d, n, v2, v1);
        elseif ty2 == 'f'
            vr = add_df(d, n, v1, v2);
        end
        
    elseif ty1 == 'f'
        if ty2 == 's'
            vr = add_sf(d, n, v2, v1);
        elseif ty2 == 'd'
            vr = add_df(d, n, v2, v1);
        end
    end
    
end


function vr = add_sd(d, n, v1, v2)

if n == 1 || d == 1
    vr = v1 + v2;
else
    vr = bsxfun(@plus, v1, v2);
end

function vr = add_sf(d, n, v1, v2)

if d == 1
    if n == 1
        vr = v1 + v2;
    else
        vr = reshape(v1, [1, 1, n]) + v2;
    end
    
else
    if n == 1
        vr = adddiag(v2, v1);
    else
        r1 = adddiag(v2(:,:,1), v1(1));
        vr = zeros(d, d, n, class(r1));
        for i = 1 : n
            vr(:,:,i) = adddiag(v2(:,:,i), v1(i));
        end
    end
end


function vr = add_df(d, n, v1, v2)

if n == 1
    if d == 1
        vr = v1 + v2;
    else
        vr = adddiag(v2, v1);
    end
else
    if d == 1
        vr = reshape(v1, [1, 1, n]) + v2;
    else
        r1 = adddiag(v2(:,:,1), v1(:,1));
        vr = zeros(d, d, n, class(r1));
        for i = 1 : n
            vr(:,:,i) = adddiag(v2(:,:,i), v1(:,i));
        end
    end
end
        
    

function vr = calc_dn(d, n, ty1, ty2, v1, v2)
% pre-conditions: 
%   S1.n == 1
%   S2.n == n > 1

if d == 1
    vr = v1 + v2;
    
    if ty1 == 'f' && (ty2 == 's' || ty2 == 'd')
        vr = reshape(vr, [1, 1, n]);
    end

else
    if ty1 == 's' || ty1 == 'd'
        
        if ty2 == 's' || ty2 == 'd'
            if ty1 == 's'
                vr = v1 + v2;
            else % d > 1 && ty1 == 'd'
                vr = bsxfun(@plus, v1, v2);
            end
            
        elseif ty2 == 'f'
            r1 = adddiag(v2(:,:,1), v1);
            vr = zeros(d, d, n, class(r1));
            vr(:,:,1) = r1;
            for i = 2 : n
                vr(:,:,i) = adddiag(v2(:,:,i), v1);
            end
        end
        
    elseif ty1 == 'f'
        if ty2 == 's' || ty2 == 'd'
            r1 = adddiag(v1, v2(:,1));
            vr = zeros(d, d, n, class(r1));
            vr(:,:,1) = r1;
            for i = 2 : n
                vr(:,:,i) = adddiag(v1, v2(:,i));
            end
            
        elseif ty2 == 'f'
            vr = bsxfun(@plus, v1, v2);
        end
        
    end
end


    