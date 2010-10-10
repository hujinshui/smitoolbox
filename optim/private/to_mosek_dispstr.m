function s = to_mosek_dispstr(level)
% convert from level to MOSEK Display string
%

switch level
    case 0
        s = 'off';
    case 1
        s = 'final';
    case 2
        s = 'iter';
    otherwise
        error('Invalid print level for mosek solver.');
end

        