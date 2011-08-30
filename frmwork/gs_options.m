function opt = gs_options(opt, varargin)
% Verify or set Gibbs sampling options
%
%   opt = gs_options(opt);
%       verify the validity of an gibbs sampling option
%
%   opt = gs_options([]);
%       returns the default gibbs sampling option struct.
%
%   opt = gs_options([], 'name1', val1, 'name2', val2, ...);
%       constructs an option struct with a series of name/value pairs.
%
%   opt = gs_options(opt, 'name1', val1, 'name2', val2, ...);
%       updates an option struct with a series of name/value pairs
%
%   Available options
%       - 'burnin':         the number of cycles for burn-in
%
%       - 'nsamples':       the number of samples to acquire 
%
%       - 'cps':            the number of cycles per sample
%
%       - 'nrepeats':       the number of repeating runs for each
%                           initial config.
%
%       - 'output_vars':    the cell array of variables to be output
%                           to sample, or an empty array to indicate
%                           output all.
%
%       - 'check':          whether to do variable checking in run-time.                           
%
%       - 'display:         the level of information displaying
%                           (can be either a string or a level number)
%                               0 | 'off':      no display
%                               1 | 'final':    display at the end
%                               2 | 'stage':    display per stage
%                               3 | 'sample':   display per sample collection
%                               4 | 'iter':     display per iteration
%

% Created by Dahua Lin, on Aug 24, 2011
%

%% main

if isempty(opt)
    opt = make_default_opt();
    
else
    if ~(isstruct(opt) && isscalar(opt) && ...
            isfield(opt, 'tag') && isequal(opt.tag, 'gs_option'))
        error('gs_options:invalidarg', ...
            'The input option struct is invalid.');
    end
end

if isempty(varargin)
    return;
end

names = varargin(1:2:end);
vals = varargin(2:2:end);

if ~(numel(names) == numel(vals) && iscellstr(names))
    error('gs_options:invalidarg', ...
        'The name value list is invalid.');
end


is_int = @(x) ...
    isscalar(x) && isnumeric(x) && isreal(x) && x == fix(x);

for i = 1 : numel(names)
    
    cn = names{i};
    lcn = lower(cn);
    cv = vals{i};
    
    switch lcn
        case {'burnin', 'nsamples', 'cps', 'nrepeats'}
            if ~(is_int(cv) && cv >= 1)
                error('gs_options:invalidarg', ...
                    'The value of option %s must be a positive integer scalar.', lcn);
            end
            opt.(lcn) = cv;
            
        case 'output_vars'
            if isempty(cv)
                opt.output_vars = [];
                
            else
                if ~iscellstr(cv)                            
                    error('gs_options:invalidarg', ...
                        'output_vars should be either empty or a cell array of strings.');
                end
                opt.output_vars = cv;
            end
            
        case 'check'
            if (isnumeric(cv) || islogical(cv)) && isscalar(cv)
                opt.check = logical(cv);
            else
                error('gs_options:invalidarg', ...
                    'check should be a logical scalar.');
            end
            
        case 'display'
            if isnumeric(cv)
                if ~(is_int(cv) && cv >= 0 && cv <= 4)
                    error('gs_options:invalidarg', ...
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
                    case 'sample'
                        dv = 3;
                    case 'iter'
                        dv = 4;
                    otherwise
                        error('gs_options:invalidarg', ...
                            'Invalid display level name ''%s''.', cv);
                end
                opt.display = dv;
            else
                error('gs_options:invalidarg', ...
                    'The value for display option is invalid.');
            end
            
        otherwise
            error('gs_options:invalidarg', ...
                'Unknown option name %s', cn);
        
    end
    
end

 
        
function opt = make_default_opt()

opt.tag = 'gs_option';
opt.burnin = 500;
opt.nsamples = 100;
opt.cps = 50;
opt.nrepeats = 1;
opt.output_vars = [];
opt.check = false;
opt.display = 0;



