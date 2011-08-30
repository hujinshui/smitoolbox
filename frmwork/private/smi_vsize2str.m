function s = smi_vsize2str(vsize)
% convert variable size spec to a string

nd = numel(vsize);
if nd == 1
    s = sprintf('vector [%s]', sn2s(vsize));
    
elseif nd == 2
    s = sprintf('matrix [%s x %s]', sn2s(vsize(1)), sn2s(vsize(2)));
    
elseif nd == 3
    s = sprintf('cube [%d x %s x %s]', sn2s(vsize(1)), sn2s(vsize(2)), sn2s(vsize(3)));
    
else
    s = 'multi-dim array (ndims >= 4)';
    
end


function s = sn2s(n)

if n >= 0    
    s = int2str(n);
else
    s = '*';
end

