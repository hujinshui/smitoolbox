function [l, u] = check_bnds(d, l, u, fname)
% Check the validity of bounds
%
%   [l, u] = check_bnds(d, l, u, fname);
%

if isempty(l)
    l = [];
else
    if ~(isnumeric(l) && isreal(l) && (isscalar(l) || isequal(size(l), [d 1])))
        error([fname ':invalidconstr'], ...
            'l should be a real column vector of size d x 1');
    end
    if ~isa(l, 'double')
        l = double(l);
    end
    if isscalar(l)
        l = constmat(d, 1, l);
    end
end

if isempty(u)
    u = [];
else
    if ~(isnumeric(u) && isreal(u) && (isscalar(u) || isequal(size(u), [d 1])))
        error([fname ':invalidconstr'], ...
            'u should be a real column vector of size d x 1');
    end
    if ~isa(u, 'double')
        u = double(u);
    end
    if isscalar(u)
        u = constmat(d, 1, u);
    end
end

