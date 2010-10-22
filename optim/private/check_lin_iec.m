function [A, bu, bl] = check_lin_iec(d, A, b, fname)
% A function to check the validity of linear inequality constraints
%
%   check_lin_iec(d, A, b, fname);
%

if isempty(A)
    if ~isempty(b)
        error([fname ':invalidconstr'], ...
            'b must be empty when A is.');
    end
    
    Aeq = [];
    bu = [];
    bl = [];
else
    if ~(isnumeric(A) && ndims(A) == 2 && isreal(A) && size(A,2) == d)
        error([fname ':invalidconstr'], ...
            'A must be a real matrix with size(A,2) == d.');               
    end
    if ~isa(A, 'double')
        A = double(A);
    end
    
    if isempty(b)
        error([fname ':invalidconstr'], ...
            'b must be non-empty when A is.');
        
    elseif iscell(b) && numel(b) == 2
        
        bl = b{1};
        bu = b{2};
        
        if isempty(bl) && isempty(bu)
            error([fname ':invalidconstr'], ...
                'bl and bu cannot be both empty.');
        end
        
    elseif isnumeric(b)
        
        bl = [];
        bu = b;
        
    else
        error([fname ':invalidconstr'], ...
            'b is in an invalid form.');
    end
    
    m = size(A, 1);
    if ~isempty(bl)
        if ~(isnumeric(bl) && isequal(size(bl), [m 1]))
            error([fname ':invalidconstr'], ...
                'bl should be an m x 1 numeric vector.');                        
        end
        if ~isa(bl, 'double')                
            bl = double(bl);
        end
    end
    
    if ~isempty(bu)
        if ~(isnumeric(bu) && isequal(size(bu), [m 1]))
            error([fname ':invalidconstr'], ...
                'bu should be an m x 1 numeric vector.');                        
        end
        if ~isa(bu, 'double')                
            bu = double(bu);
        end
    end
end

       