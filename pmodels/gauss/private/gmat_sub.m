function R = gmat_sub(cf, C, si)
% Get a subset of the matrice(s)
%

switch cf
    case 's'
        R = C(1, si);
    case 'd'
        R = C(:, si);
    case 'f'
        R = C(:,:,si);
end