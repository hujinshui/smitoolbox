function v = get_display_level(funname, s)
% convert the display level from string to value
%
%   v = get_display_level(funname, s);
%

switch s
    case 'off'
        v = 0;
    case 'notify'
        v = 1;
    case 'final'
        v = 2;
    case 'iter'
        v = 3;
    otherwise
        error([funname ':invalidopt'], 'Invalid display level %s', s);
end
