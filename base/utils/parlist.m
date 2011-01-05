function S = parlist(varargin)
% Parses a name/value pair list
%
%   S = parlist('name1', value1, 'name2', value2, ...)
%       parses a name/value list and return the results as a struct S.
%
%       S is a struct with S.name1 = value1, S.name2 = value2, ...
%
%       Note that the value names should be valid variable name that can
%       serve as a field name in S.
%
%   S = parlist(S0, 'name1', value1, 'name2', value2, ...)
%       updates the struct S0 with the values provided in the ensuing
%       name/value list.
%
%       Each name in the name pair list should be a field of S0.
%
%   Remarks
%   -------
%       - If there are repeated names in the list, then the latter value
%         corresponding to the same name will override the former one.
%

%   History
%   -------
%       - Created by Dahua Lin, on Nov 10, 2010
%       - Modified by Dahua Lin, on Jan 4, 2011
%

%% main

if nargin == 0
    S = [];
    
else
    if ischar(varargin{1})
        S = make_s([], varargin);        
    elseif isstruct(varargin{1})
        S = make_s(varargin{1}, varargin(2:end));        
    end
    
end

%% core 

function S = make_s(S, pairs)
% the core function to make the struct

if ~isempty(pairs)

    names = pairs(1:2:end);
    values = pairs(2:2:end);

    n = numel(names);
    if ~(n == numel(values) && iscellstr(names))
        error('parlist:invalidarg', 'the name/value pair list is invalid.');
    end
    
    for i = 1 : n
        cname = names{i};
        if ~isfield(S, cname)
            error('parlist:invalidarg', 'The option name %s is not supported', cname);
        end
        S.(cname) = values{i};
    end
end

