function tf = smi_is_valid_typename(name)
% Tests whether a name is a valid typename

if isvarname(name)
    tf = exist(name, 'class') || ...
        strcmp(name, 'numeric') || strcmp(name, 'float'); 
else
    tf = false;
end
