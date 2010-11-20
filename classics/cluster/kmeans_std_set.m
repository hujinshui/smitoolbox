function S = kmeans_std_set(S0, varargin)
% Set options for kmeans_std function.
%
%   S = kmeans_std_set;
%       returns the default option struct.
%
%   S = kmeans_std_set(S0, 'name1', value1, 'name2', value2, ...);
%       Verifies input options and constructs an option struct for 
%       kmeans_ex_set.
%
%       Here, S0 is the original option struct to be updated. 
%       It can be empty, when constructing based on default options.
%
%       Please refer to the help of kmeans_ex for details.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%

%% main

% default options

if nargin < 1 || isempty(S0)    
    S = struct( ...
        'tag', 'kme_std_checked', ...
        'maxiter', 100, ...
        'tolc', 0, ...
        'display', 'off', ...
        'displevel', 0, ...
        'ucwarn', false, ...
        'init', 'km++', ...
        'initfunc', @kmpick_pp, ...
        'onnil', 'repick++', ...
        'onnilop', 1, ...
        'rpickfunc', @kmpick_pp, ...
        'distfunc', @kmd_sqL2, ...
        'rstream', []);    
else    
    if ~(isstruct(S0) && numel(S0) == 1 && isfield(S0, 'tag') && ...
            strcmp(S0, 'kme_std_checked'))
        error('kmeans_ex_set:invalidarg', 'S0 is invalid.');
    end    
    S = S0;
end
    
% update struct

if isempty(varargin)
    return;
end

onames = varargin(1:2:end);
ovals = varargin(2:2:end);

n = length(onames);
if ~(iscellstr(onames) && n == length(ovals))
    error('kmeans_std_set:invalidarg', 'The name/value pair list is invalid.');
end

update_display = false;
update_init = false;
update_on_nil = false;

for i = 1 : n
   
    name = lower(onames{i});
    v = ovals{i};
    
    switch name
        
        case {'maxiter', 'tolc'}
            if ~(isnumeric(v) && isscalar(v) && v > 0)
                opterr('%s must be a positive scalar.', name);
            end
                        
        case 'display'
            if ~(ischar(v) && any(strcmpi(v, {'off', 'iter', 'final'})))
                opterr('display must be either of ''off'', ''iter'', or ''final''.');
            end
            v = lower(v);
            update_display = ~strcmp(v, S.display);                
            
        case 'ucwarn'
            if ~((islogical(v) || isnumeric(v)) && isscalar(v))
                opterr('uc_warn must be a logical or numeric scalar (0 or 1).');
            end
            v = logical(v);            
            
        case 'init'            
            if ~(ischar(v) && any(strcmpi(v, {'random', 'km++', 'mcinit'})))
                opterr('The value of init option is invalid.');
            end
            v = lower(v);
            update_init = ~strcmp(v, S.init);
            
        case 'onnil'
            if ~(ischar(v) && any(strcmpi(v, {'repick', 'repick++', 'mcpick', 'error', 'keep'})))
                opterr('The value of on_nil option is invalid.');
            end
            v = lower(v);
            update_on_nil = ~strcmp(v, S.on_nil);
            
        case 'distfunc'
            if ischar(v)                
                if ~any(strcmpi(v, {'sqL2', 'L1'}))
                    opterr('Unknown dist_func name %s', v);
                end
                v = lower(v);
                if strcmpi(v, 'sqL2')
                    v = @kmd_sqL2;
                elseif strcmpi(v, 'L1')
                    v = @kmd_L1;
                end
                
            else
                if ~isa(v, 'function_handle')
                    opterr('The value of dist_func is invalid.');
                end
            end
            
        case 'rstream'
            if ~(isempty(v) || isa(v, 'RandStream'))
                opterr('rstream should be either empty or a RandStream object.');
            end
            if isempty(v)
                v = [];
            end
            
        otherwise
            error('kmeans_ex_set:invalidarg', 'Invalid option name %s', name);
    end
    
    S.(name) = v;
end


% process options

if update_display
    switch S.display
        case 'off'
            S.displevel = 0;
        case 'final'
            S.displevel = 1;
        case 'iter'
            S.displevel = 2;
    end
end

if update_init
    switch S.init
        case 'km++'
            S.initfunc = @kmpick_pp;
        case 'random'
            S.initfunc = @kmpick_rand;
        case 'mcinit'
            S.initfunc = @kmpick_mc;
    end
end

if update_on_nil
    switch S.onnil
        case 'repick++'
            S.rpickfunc = @kmpick_pp;
            S.onnilop = 1;
        case 'repick'
            S.rpickfunc = @kmpick_rand;
            S.onnilop = 1;
        case 'mcpick'
            S.rpickfunc = @kmpick_mc;
            S.onnilop = 1;
        case 'error'
            S.rpickfunc = [];
            S.onnilop = -1;
        case 'keep'
            S.rpickfunc = [];
            S.onnilop = 0;
    end
end


%% sub function

function opterr(msg, varargin)

error('kmeans_std_set:invalidopt', msg, varargin{:});


