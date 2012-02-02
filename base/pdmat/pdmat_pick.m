function R = pdmat_pick(S, i)
% Get a subset from a pack of positive definite matrices
%
%   R = pdmat_pick(S, i);
%       returns a pdmat struct that is comprised of the i-th matrix
%       in S of the same form.
%
%       i can also be a index vector or logical vector, in which case,
%       R may contain multiple matrices.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% main

switch S.ty
    case 's'
        v = S.v(1, i);
        n = size(v, 2);
    case 'd'       
        v = S.v(:, i);
        n = size(v, 2);
    case 'f'
        v = S.v(:,:,i);
        n = size(v, 3);
end

R.tag = S.tag;
R.ty = S.ty;
R.d = S.d;
R.n = n;
R.v = v;
        