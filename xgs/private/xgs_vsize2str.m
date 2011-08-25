function s = xgs_vsize2str(vsize)
% convert variable size spec to a string

nd = numel(vsize);
if nd == 1
    s = sprintf('vector [%d]', vsize);
    
elseif nd == 2
    s = sprintf('matrix [%d x %d]', vsize(1), vsize(2));
    
elseif nd == 3
    s = sprintf('cube [%d x %d x %d]', vsize(1), vsize(2), vsize(3));
    
else
    s = 'multi-dim array (ndims >= 4)';
    
end