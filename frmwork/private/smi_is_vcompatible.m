function tf = smi_is_vcompatible(vs, vd)
% test whether the variable described in vs can be assigned to one as vd
%

% check type

tf = false;

if strcmp(vs.type, vd.type)
    tf = true;
    
elseif strcmp(vd.type, 'any')
    tf = true;
    
else
    switch vs.type
        case {'double', 'single'}
            s_isfloat = 1;
            s_isnum = 1;
        case {'uint8', 'int8', 'uint16', 'int16'}
            s_isfloat = 0;
            s_isnum = 1;
        case {'uint32', 'int32', 'uint64', 'int64'}
            s_isfloat = 0;
            s_isnum = 1;
        otherwise
            s_isfloat = 0;
            s_isnum = 0;
    end
    
    if s_isfloat && strcmp(vd.type, 'float')
        tf = true;
    elseif s_isnum && strcmp(vd.type, 'numeric')
        tf = true;
    end
end

if ~tf
    return;
end        

% check size

nd = numel(vs.size);
if nd ~= numel(vd.size)
    tf = false;
    return;
end

tf = all(vd.size == vs.size | vd.size == -1);
