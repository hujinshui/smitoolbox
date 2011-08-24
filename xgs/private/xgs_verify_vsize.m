function vsize = xgs_verify_vsize(vsize)
% Check whether the input is a valid variable size spec
%
%  If valid, returns a processed vsize;
%  otherwise, returns empty.
%

if ~(isnumeric(vsize) && ~isempty(vsize) && ndims(vsize) == 2 && size(vsize,1) == 1)
    vsize = [];
    return;
end

if ~(all(vsize >= 0 & vsize == fix(vsize)))
    vsize = [];
    return;
end

if all(vsize == 1)
    vsize = 1;
else
    nd = find(vsize ~= 1, 1, 'last');
    vsize = vsize(1:nd);
    if ~isa(vsize, 'double')
        vsize = double(vsize);
    end
end
