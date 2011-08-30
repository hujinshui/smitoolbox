function M = smi_make_map(L)
% Makes name->id map from a list of named items

M = [];
n = numel(L);

if isstruct(L)
    for i = 1 : n
        M.(L(i).name) = i;
    end
elseif iscell(L)        
    for i = 1 : n
        M.(L{i}.name) = i;
    end
end
