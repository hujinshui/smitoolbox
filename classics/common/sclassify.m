function L = sclassify(V, op)
% Do standard classification (assign the label with min or max value)
%
%   L = sclassify(V, 'max');
%       Assigns each column a label whose value is the maximum.
%
%       Let V be an m x n matrix, then L is a vector of size 1 x n.
%       Such that V(L(i), i) is the maximum of the values in column
%       V(:, i).
%
%   L = sclassify(V, 'min');
%       Assigns each column a label whose values is the minimum.
%

% Created by Dahua Lin, on Jun 6, 2010
%

%% verify input

if ~(isnumeric(V) && ndims(V) == 2)
    error('sclassify:invalidarg', ...
        'V should be a numeric matrix.');
end

if ~ischar(op)
    error('sclassify:invalidarg', 'The second arg should be a char string.');
end

if strcmp(op, 'max')
    is_max = true;
elseif strcmp(op, 'min')
    is_max = false;
else
    error('sclassify:invalidarg', 'The second arg should be either ''max'' or ''min''.');
end
   
%% main

if is_max    
    [dummy, L] = max(V, [], 1); %#ok<ASGLU>
else
    [dummy, L] = min(V, [], 1); %#ok<ASGLU>
end

    
    
