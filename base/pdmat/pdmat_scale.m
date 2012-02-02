function Sr = pdmat_scale(S, c)
% Compute a scalar-product of a positive definite matrix
%
%   Sr = pdmat_scale(S, c);
%       returns Sr, whose matrix is the scalar product of c and
%       that in S.
%
%       c can be a scalar, or a row vector. 
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

if ~(isfloat(c) && isreal(c))
    error('pdmat_scale:invalidarg', ...
        'c should be of real floating-point type.');
end

if ~(isscalar(c) || (ndims(c) == 2 && size(c,1) == 1))
    error('pdmat_scale:invalidarg', ...
        'c should be either a scalar or a row vector.');
end

if S.n > 1
    K = size(c, 2);
    if ~(K == 1 || K == S.n)
        error('pdmat_scale:invalidarg', ...
            'When S.n > 1, length(c) should either 1 or S.n.');
    end
end


%% main

if isscalar(c)
    Sr = S;
    Sr.v = c * S.v;
    
else % size(c, 2) = K > 1
    
    Sr = S;    
    ty = S.ty;
    d = S.d;
        
    Sr.n = size(c, 2);
    
    switch ty
        case 's'
            Sr.v = S.v .* c;
        case 'd'
            if d == 1
                Sr.v = S.v .* c;
            else
                Sr.v = bsxfun(@times, S.v, c);
            end
        case 'f'
            c3 = reshape(c, [1, 1, Sr.n]);
            if d == 1
                Sr.v = S.v .* c3;
            else
                Sr.v = bsxfun(@times, S.v, c3);
            end
    end   
    
end

