function tf = xgs_is_vcompatible(vs, vd)
% test whether the variable described in vs can be assigned to one as vd
%

% check type

if ~(isequal(vs.type, vd.type) || isequal(vd.type, 'any'))
    tf = false;
    return;
end

% check size

nd = numel(vs.size);
if nd ~= numel(vd.size)
    tf = false;
    return;
end

tf = all(vd.size == vs.size | vd.size == 0);
