function [A, bl, bu] = check_lin_iec(d, A, b, funname)
% Check the validity of linear inequality constraints
%
%   [A, bl, bu] = check_lin_iec(d, A, b, funname);
%

if ~(isfloat(A) && ndims(A) == 2)
    error([funname ':invalidconstr'], ...
        'The inequality coefficient matrix A should be a numeric matrix.');
end

[m, n] = size(A);
if n ~= d
    error([funname ':invalidconstr'], ...
        'The requirement size(A, 2) == d is not satisfied.');
end

if iscell(b) && numel(b) == 2
    bl = b{1};
    bu = b{2};            
    if ~(isempty(bl) || (isfloat(bl) && ndims(bl) == 2 && size(bl, 2) == 1))
        error([funname ':invalidconstr'], ...
            'The left-bound of inequality bl should be a numeric column vector.');
    end
    if ~(isempty(bu) || (isfloat(bu) && ndims(bu) == 2 && size(bu, 2) == 1))
        error([funname ':invalidconstr'], ...
            'The right-bound of inequality bu should be a numeric column vector.');
    end
    
    if ~((isempty(bl) || size(bl,1) == m)  && (isempty(bu) || size(bu,1) == m))
        error([funname ':invalidconstr'], ...
            'The size of the bounds of inequality is inconsistent with A.');
    end
    
elseif isnumeric(b)
    if ~(isfloat(b) && ndims(b) == 2 && size(b, 2) == 1)
        error([funname ':invalidconstr'], ...
            'The right-bound of inequality b should be a numeric column vector.');
    end
    if size(b,1) ~= m
        error([funname ':invalidconstr'], ...
            'The size of the bounds of inequality is inconsistent with A.');
    end
    
    bl = [];
    bu = b;
end

    
