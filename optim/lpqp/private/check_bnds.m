function [l, u] = check_bnds(d, l, u, funname)
% Check the validity of range bounds
%
%   [l, u] = check_bnds(d, l, u, funcname);
%

if ~isempty(l)
    if ~(isfloat(l) && (isscalar(l) || isequal(size(l), [d 1])))
        error([funname ':invalidconstr'], ...
            'The lower bound should be either a scalar of a d x 1 column.');
    end
    
    if isscalar(l) 
        l = constmat(d, 1, l);
    end
else
    l = [];
end

if ~isempty(u)
    if ~(isfloat(u) && (isscalar(u) || isequal(size(u), [d 1])))
        error([funname ':invalidconstr'], ...
            'The upper bound should be either a scalar of a d x 1 column.');
    end
    
    if isscalar(u)
        u = constmat(d, 1, u);
    end
else
    u = [];
end
   
