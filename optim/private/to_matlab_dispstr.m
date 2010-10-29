function s = to_matlab_dispstr(level)
% convert from level to MATLAB Display string
%

switch level
    case 0
        s = 'off';
    case 1
        s = 'notify';
    case 2
        s = 'final';
    case 3
        s = 'iter';
    otherwise
        error('Invalid print level for matlab solver.');
end

        