function opt = mcmc_options(opt, varargin)
% Verify or set MCMC sampling options for smi_mcmc
%
%   opt = mcmc_options();
%   opt = mcmc_options([]);
%       returns the default MCMC sampling option struct.
%
%   opt = mcmc_options(opt);
%       verify the validity of an MCMC sampling option
%
%   opt = mcmc_options([], 'name1', val1, 'name2', val2, ...);
%       constructs an option struct with a series of name/value pairs.
%
%   opt = mcmc_options(opt, 'name1', val1, 'name2', val2, ...);
%       updates an option struct with a series of name/value pairs
%
%   Available options
%       - 'burnin':         the number of iterations for burn-in
%
%       - 'nsamples':       the number of samples to acquire 
%
%       - 'ips':            the number of iterations per sample         
%
%       - 'display:         the level of information displaying
%                           (can be either a string or a level number)
%                               0 | 'off':      no display
%                               1 | 'final':    display at the end
%                               2 | 'stage':    display per stage
%                               3 | 'sample':   display per sample collection
%                               4 | 'iter':     display per iteration
%

% Created by Dahua Lin, on Sep 4, 2011
%

%% main

if nargin < 1 || isempty(opt)
    opt = make_default_opt();
    
else
    if ~(isstruct(opt) && isscalar(opt) && ...
            isfield(opt, 'tag') && isequal(opt.tag, 'mcmc_options'))
        error('mcmc_options:invalidarg', ...
            'The input option struct is invalid.');
    end
end

if isempty(varargin)
    return;
end

names = varargin(1:2:end);
vals = varargin(2:2:end);

if ~(numel(names) == numel(vals) && iscellstr(names))
    error('mcmc_options:invalidarg', ...
        'The name value list is invalid.');
end


is_int = @(x) ...
    isscalar(x) && isnumeric(x) && isreal(x) && x == fix(x);

for i = 1 : numel(names)
    
    cn = names{i};
    lcn = lower(cn);
    cv = vals{i};
    
    switch lcn
        case {'burnin', 'nsamples', 'ips'}
            if ~(is_int(cv) && cv >= 1)
                error('mcmc_options:invalidarg', ...
                    'The value of option %s must be a positive integer scalar.', lcn);
            end
            opt.(lcn) = cv;
                        
        case 'display'
            if isnumeric(cv)
                if ~(is_int(cv) && cv >= 0 && cv <= 4)
                    error('mcmc_options:invalidarg', ...
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
                        error('mcmc_options:invalidarg', ...
                            'Invalid display level name ''%s''.', cv);
                end
                opt.display = dv;
            else
                error('mcmc_options:invalidarg', ...
                    'The value for display option is invalid.');
            end
            
        otherwise
            error('mcmc_options:invalidarg', ...
                'Unknown option name %s', cn);
        
    end
    
end

 
        
function opt = make_default_opt()

opt.tag = 'mcmc_options';
opt.burnin = 500;
opt.nsamples = 100;
opt.ips = 50;
opt.nrepeats = 1;
opt.display = 0;



