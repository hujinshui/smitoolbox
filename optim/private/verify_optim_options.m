function verify_optim_options(funname, options, varargin)
% verify specified optimization options
%
%   verify_optim_options(funname, options, <name1>, <name2>, ...);
%

if ~isstruct(options) && numel(options) == 1
    error(funname + ':invalidopt', 'options should be a struct.');
end

n = length(varargin);
for i = 1 : n
    
    oname = varargin{i};
    
    if ~isfield(options, oname)
        error([funname ':invalidopt'], ...
            'The options lacks the field %s', oname);
    end
    
    v = options.(oname);
    
    switch oname
        case 'MaxIter'
            if ~(isnumeric(v) && isscalar(v) && v >= 1 && v == fix(v))
                error([funname ':invalidopt'], ...
                    '%s should be a positive integer', oname);                
            end
            
        case {'TolX', 'TolFun'}
            if ~(isfloat(v) && isreal(v) && isscalar(v) && v > 0)
                error([funname ':invalidopt'], ...
                    '%s should be a positive real value', oname);
            end
            
        case 'Display'
            if ~ischar(v)
                error([funname ':invalidopt'], ...
                    '%s should be a char string', oname);
            end
    end    
    
end

    
    