function opt = xgs_options(opt, varargin)
% Verify or set XGS options
%
%   opt = xgs_options(opt);
%       verify the validity of an XGS option
%
%   opt = xgs_options([]);
%       returns the default xgs option struct.
%
%   opt = xgs_options([], 'name1', val1, 'name2', val2, ...);
%       constructs an xgs option struct with a series of name/value pairs.
%
%   opt = xgs_options(opt, 'name1', val1, 'name2', val2, ...);
%       updates an xgs option struct with a series of name/value pairs
%
%   Available options
%       - 'burnin':         the number of cycles for burn-in
%
%       - 'nsamples':       the number of samples to acquire 
%
%       - 'cps':            the number of cycles per sample
%
%       - 'output_vars':    the cell array of variables to be output
%                           to sample, or an empty array to indicate
%                           output all.
%
%       - 'display:         the level of information displaying
%                           (can be either a string or a level number)
%                               0 | 'off':      no display
%                               1 | 'final':    display at the end
%                               2 | 'stage':    display per stage
%                               3 | 'cycle':    display per cycle
%                               4 | 'step':     display per step
%

% Created by Dahua Lin, on Aug 24, 2011
%

%% main

if isempty(opt)
    opt = make_default_opt();
    
else
    if ~(isstruct(opt) && isscalar(opt) && ...
            isfield(opt, 'tag') && isequal(opt.tag, 'xgs_option'))
        error('xgs_options:invalidarg', ...
            'The input option struct is invalid.');
    end
end

if isempty(varargin)
    return;
end

names = varargin(1:2:end);
vals = varargin(2:2:end);

if ~(numel(names) == numel(vals) && iscellstr(names))
    error('xgs_options:invalidarg', ...
        'The name value list is invalid.');
end


is_int = @(x) ...
    isscalar(x) && isnumeric(x) && isreal(x) && x == fix(x);

for i = 1 : numel(names)
    
    cn = names{i};
    lcn = lower(cn);
    cv = vals{i};
    
    switch lcn
        case {'burnin', 'nsamples', 'cps'}
            if ~(is_int(cv) && cv >= 1)
                error('xgs_options:invalidarg', ...
                    'The value of option %s must be a positive integer scalar.', lcn);
            end
            opt.(lcn) = cv;
            
        case 'output_vars'
            if isempty(cv)
                opt.output_vars = [];
                
            else
                if ~iscellstr(cv)                            
                    error('xgs_options:invalidarg', ...
                        'output_vars should be either empty or a cell array of strings.');
                end
                opt.output_vars = cv;
            end
            
        case 'display'
            if isnumeric(cv)
                if ~(is_int(cv) && cv >= 0 && cv <= 4)
                    error('xgs_options:invalidarg', ...
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
                    case 'cycle'
                        dv = 3;
                    case 'step'
                        dv = 4;
                end
                opt.display = dv;
            else
                error('xgs_options:invalidarg', ...
                    'The value for display option is invalid.');
            end
            
        otherwise
            error('xgs_options:invalidarg', ...
                'Unknown XGS option name %s', cn);
        
    end
    
end

 
        
function opt = make_default_opt()

opt.tag = 'xgs_option';
opt.burnin = 500;
opt.nsamples = 100;
opt.cps = 50;
opt.output_vars = [];
opt.display = 0;



