function tf = smi_verify_var(v, vspec)
% Tests the validity of a variable against a variable specification

tf = true;

% verify type

if ~isa(v, vspec.type)
    tf = false;
    return;
end

% verify size

vsiz = vspec.size;
nd = numel(vsiz);

if nd == 1
    if ~(ndims(v) == 2 && size(v,2) == 1)
        tf = false;
        return;
    end
    
    if (vsiz >= 0 && size(v,1) ~= vsiz)
        tf = false;
        return;
    end
    
else % nd > 1
    
    if ndims(v) ~= nd
        tf = false;
        return;
    end
    
    asiz = size(v);
    
    if ~all(asiz == vsiz | vsiz < 0)
        tf = false;
        return;
    end
end

