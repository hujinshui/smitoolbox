function Sr = pdmat_scale(S, c)
% Compute a scalar-product of a positive definite matrix
%
%   Sr = pdmat_scale(S, c);
%       returns Sr, whose matrix is the scalar product of c and
%       that in S.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% verify input

if ~(isfloat(c) && isscalar(c) && isreal(c))
    error('pdmat_scale:invalidarg', ...
        'c should be a numeric scalar.');
end

%% main

Sr = S;
Sr.v = c * S.v;

