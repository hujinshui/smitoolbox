function v = uniqvalues(x)
% Get a list of sorted values with no repetition
%
%   v = uniqvalues(x);
%       returns a row vector v that contains all values in array x with
%       no repetition.
%
%   Remarks
%   -------
%       - It is similar to the unique function shiped with MATLAB. 
%         However, it is implemented with C++ mex, and specifically
%         designed for numeric arrays. Hence, it is much more efficient.
%

% Created by Dahua Lin, on June 6, 2010
%


%% verify input

if ~isnumeric(x) || issparse(x)
    error('uniqvalues:invalidarg', ...
        'x should be a non-sparse numeric array.');
end

%% main

sx = sort(reshape(x, 1, numel(x)));
p = valueseg(sx);

v = sx(p);

