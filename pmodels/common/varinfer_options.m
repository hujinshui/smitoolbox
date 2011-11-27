function opt = varinfer_options(opt, varargin)
% Verify or set variational inference control options for smi_varinfer
%
%   opt = varinfer_options();
%   opt = varinfer_options([]);
%       returns the default variational inference option struct.
%
%   opt = varinfer_options(opt);
%       verify the validity of an variational inference option
%
%   opt = varinfer_options([], 'name1', val1, 'name2', val2, ...);
%       constructs an option struct with a series of name/value pairs.
%
%   opt = varinfer_options(opt, 'name1', val1, 'name2', val2, ...);
%       updates an option struct with a series of name/value pairs
%
%   Available options
%       - 'maxiters':       the maximum number of iterations
%
%       - 'ipe':            the number of iterations between two
%                           evaluations of the objective value.
%
%       - 'tol':            the tolerance of objective change between
%                           current and last objective values on
%                           convergence.
%
%       - 'nrepeats':       the number of repeating runs. If using
%                           random initialization, different runs
%                           may yield different results.
%
%       - 'display:         the level of information displaying
%                           (can be either a string or a level number)
%                           0 | 'off':      no display
%                           1 | 'final':    display at the end
%                           2 | 'stage':    display per stage
%                           3 | 'eval':     display per objective evaluation
%                           4 | 'iter':     display per iteration
%

% Created by Dahua Lin, on Sep 4, 2011
%

%% main

if nargin < 1 || isempty(opt)
    opt = make_default_opt();
    
else
    if ~(isstruct(opt) && isscalar(opt) && ...
            isfield(opt, 'tag') && isequal(opt.tag, 'varinfer_options'))
        error('varinfer_options:invalidarg', ...
            'The input option struct is invalid.');
    end
end

if isempty(varargin)
    return;
end

names = varargin(1:2:end);
vals = varargin(2:2:end);

if ~(numel(names) == numel(vals) && iscellstr(names))
    error('varinfer_options:invalidarg', ...
        'The name value list is invalid.');
end


is_int = @(x) ...
    isscalar(x) && isnumeric(x) && isreal(x) && x == fix(x);

for i = 1 : numel(names)
    
    cn = names{i};
    lcn = lower(cn);
    cv = vals{i};
    
    switch lcn
        case {'maxiters', 'ipe', 'nrepeats'}
            if ~(is_int(cv) && cv >= 1)
                error('varinfer_options:invalidarg', ...
                    'The value of option %s must be a positive integer scalar.', lcn);
            end
            opt.(lcn) = cv;
            
        case 'tol'
            if ~(isfloat(cv) && isscalar(cv) && cv >= 0)
                error('varinfer_options:invalidarg', ...
                    'The value of option tol must be a non-negative real number.');
            end
            opt.tol = cv;
                        
        case 'display'
            if isnumeric(cv)
                if ~(is_int(cv) && cv >= 0 && cv <= 4)
                    error('varinfer_options:invalidarg', ...
                        'The value of cps must be an integer in [0, 4].');
                end
                opt.display = cv;
                
            elseif ischar(cv)
                switch lower(cv)
                    case 'off'
                        dv = 0;
                    case 'final'
                        dv = 1;
                    case 'stage'
                        dv = 2;
                    case 'eval'
                        dv = 3;
                    case 'iter'
                        dv = 4;
                    otherwise
                        error('varinfer_options:invalidarg', ...
                            'Invalid display level name ''%s''.', cv);
                end
                opt.display = dv;
            else
                error('varinfer_options:invalidarg', ...
                    'The value for display option is invalid.');
            end
            
        otherwise
            error('varinfer_options:invalidarg', ...
                'Unknown option name %s', cn);
        
    end
    
end

 
        
function opt = make_default_opt()

opt.tag = 'varinfer_options';
opt.maxiters = 200;
opt.ipe = 1;
opt.tol = 1e-6;
opt.nrepeats = 1;
opt.display = 0;



