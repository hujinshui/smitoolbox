function [Aeq, beq] = check_lin_eqc(d, Aeq, beq, fname)
% A function to check the validity of linear equality constraints
%
%   check_lin_eqc(d, Aeq, beq, fname);
%

if isempty(Aeq)
    if ~isempty(beq)
        error([fname ':invalidconstr'], ...
            'beq must be empty when Aeq is.');
    end
    
    Aeq = [];
    beq = [];    
else
    if ~(isnumeric(Aeq) && ndims(Aeq) == 2 && isreal(Aeq) && size(Aeq,2) == d)
        error([fname ':invalidconstr'], ...
            'Aeq must be a real matrix with size(Aeq,2) == d.');               
    end
    if ~isa(Aeq, 'double')
        Aeq = double(Aeq);
    end
    
    if isempty(beq)
        error([fname ':invalidconstr'], ...
            'beq must be non-empty when A is.');
                
    elseif isnumeric(beq)
        
        m = size(Aeq, 1);
        if ~(isreal(beq) && isequal(size(beq), [m 1]))
            error([fname ':invalidconstr'], ...
                'bl should be an m_eq x 1 numeric vector.');                        
        end
        if ~isa(beq, 'double')                
            beq = double(beq);
        end
        
    else
        error([fname ':invalidconstr'], ...
            'beq is in an invalid form.');
    end
        
end