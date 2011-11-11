function M = pdmat_fullform(S, i)
% Gets the full matrix form of a matrix packed in pdmat
%
%   M = pdmat_fullform(S);
%       Gets the full matrix representation of the matrix represented
%       by a pdmat struct S.
%
%       Here, S should contain only one matrix (i.e. S.n == 1)
%
%   M = pdmat_fullform(S, i);
%       Gets the full matrix representation of the i-th matrix in S.
%

% Created by Dahua Lin, on Aug 25, 2010
%

%% main

% verify input and get v (of the selected part)

if nargin < 2
    if S.n ~= 1
        error('pdmat_fullform:invalidarg', ...
            'S should contain only one matrix without i given.');
    end
    
    v = S.v;
    
else
    if ~(isnumeric(i) && isscalar(i) && i == fix(i) && i >= 1 && i <= S.n)
        error('pdmat_fullform:invalidarg', ...
            'i should be an integer scalar in [1, S.n].');
    end
    
    switch S.ty
        case 's'
            v = S.v(i);
        case 'd'
            v = S.v(:, i);
        case 'f'
            v = S.v(:,:,i);
    end
end

% make the full matrix

switch S.ty
    case 's'
        M = diag(v * ones(S.d, 1));
    case 'd'
        M = diag(v);
    case 'f'
        M = v;
end


