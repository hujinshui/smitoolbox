function r = xgs_from_repr(repr, I)

n = numel(I);
r = cell(1, n);
for i = 1 : n
    if I(i) > 0
        r{i} = repr{I(i)};
    end
end
