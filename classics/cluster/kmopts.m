function S = kmopts(varargin)
%KMOPTS Options for K-means and K-medoid
%
%   S = KMOPTS;
%
%       returns the default option struct.
%
%   S = KMOPTS('name1', value1, 'name2', value2, ...);
%
%       constructs an option struct, modifying part of the option values
%       using the input values.
%
%   S = KMOPTS(S0, 'name1', value1, 'name2', value2, ...);
%
%       constructs a new option struct, updating the values in an old
%       option struct given by S0.
%
%       Please refer to the help of kmeans_ex for details.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 27, 2010
%       - Modified by Dahua Lin, on Mar 27, 2012
%

%% parse inputs

if nargin == 0
    S = default_opts();
    return;
    
elseif ischar(varargin{1})
    S = default_opts();
    nvlist = varargin;
    
elseif isstruct(varargin{1})
    S = varargin{1};
    if ~(isscalar(S) && isfield(S, 'tag') && strcmp(S.tag, 'kmopts'))
        error('kmopts:invalidarg', 'The struct S0 is invalid.');
    end
    if numel(varargin) == 1
        return;
    end    
    nvlist = varargin(2:end);
    
else
    error('kmopts:invalidarg', 'Invalid input arguments.');
end
    


%% main
    
% update struct

onames = nvlist(1:2:end);
ovals = nvlist(2:2:end);

n = length(onames);
if ~(iscellstr(onames) && n == length(ovals))
    error('kmopts:invalidarg', ...
        'The name/value pair list is invalid.');
end


for i = 1 : n
   
    name = lower(onames{i});
    v = ovals{i};
    
    switch name
        
        case 'maxiter'
            if ~(isnumeric(v) && isscalar(v) && v > 0)
                opterr('%s must be a positive scalar.', name);
            end
            
        case {'tolc', 'tolfun'}
            if ~(isnumeric(v) && isscalar(v) && isreal(v) && v >= 0)
                opterr('%s must be a non-negative real scalar.', name);
            end                        
                        
        case 'display'
            if ~ischar(v)
                opterr('display must be a char string.');
            end
            v = lower(v);
                         
            switch v
                case 'off'
                    S.displevel = 0;
                case 'final'
                    S.displevel = 1;
                case 'iter'
                    S.displevel = 2;
                otherwise
                    error('The value of display is invalid.');
            end
            
        otherwise
            error('kmopts:invalidarg', 'Invalid option name %s', name);
    end
    
    S.(name) = v;
end


%% sub function

function S = default_opts()

S = struct( ...
    'tag', 'kmopts', ...
    'maxiter', 100, ...
    'tolfun', 1e-8, ...
    'tolc', 0, ...
    'display', 'off', ...
    'displevel', 0);


function opterr(msg, varargin)

error('kmopts:invalidopt', msg, varargin{:});


