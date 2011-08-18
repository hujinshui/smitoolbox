function n = gmat_num(cf, C)
% get the number of matrices packed in an array
%

switch cf
    case 's'
        n = numel(C);
    case 'd'
        n = size(C, 2);
    case 'f'
        n = size(C, 3);
end