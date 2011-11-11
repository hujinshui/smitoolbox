function [u, varargout] = uniqvalues(x, tags)
% Get a list of sorted values with no repetition
%
%   u = uniqvalues(x);
%       returns a row vector v that contains all values in array x with
%       no repetition.
%
%   [u, I] = uniqvalues(x, 'I');
%       additionally returns an integer map, such that x(i) == u(I(i)).
%
%   [u, G] = uniqvalues(x, 'G');
%       additionally returns a cell array of grouped indices, such that
%       all(x(G{k}) == u(k)).
%
%   [u, C] = uniqvalues(x, 'C');
%       additionally returns the counts of all distinct values.
%
%   [u, ...] = uniqvalues(x, tags);
%       returns additional information according to tags.
%
%       For example, one can write the following to have the function
%       return both the grouped indices and the counts.
%
%           [u, G, C] = uniqvalues(x, 'GC');
%
%
%   Remarks
%   -------
%       - It is similar to the builtin unique function. 
%         However, it is implemented with C++ mex, and specifically
%         designed for numeric arrays. Hence, it is much more efficient.
%

%   History
%   -------
%       - Created by Dahua Lin, on June 6, 2010
%       - Modified by Dahua Lin, on Nov 11, 2010
%           - supports multiple outputs
%


%% verify input

if ~isnumeric(x) || issparse(x)
    error('uniqvalues:invalidarg', ...
        'x should be a non-sparse numeric array.');
end

xsiz = size(x);
if ~(ndims(x) == 2 && xsiz(1) == 1)
    x = reshape(x, 1, numel(x));
end


use_I = 0;
use_C = 0;
use_G = 0;

if nargin >= 2 && ~isempty(tags)
    if ~ischar(tags)
        error('uniqvalues:invalidarg', 'the tags should be a char array.');
    end
    
    if any(tags == 'C')
        use_C = 1;
    end
    
    if any(tags == 'I')
        use_C = 1;
        use_I = 1;
    end
    
    if any(tags == 'G')
        use_G = 1;
    end                
else
    tags = [];
end


%% main

% extract unique values

if use_I || use_G
    [sx, si] = sort(x);
else
    sx = sort(x);
end
p = valueseg(sx);
u = sx(p);

% additional output

if use_C
    C = [p(2:end), numel(x)+1] - p;
end

if use_I
    I = zeros(xsiz);
    I(si) = repnum(C);
end

if use_G
    K = numel(p);
    ep = [p(2:end)-1, numel(x)];
    G = cell(1, K);
    for k = 1 : K
        G{k} = si(p(k):ep(k));
    end
end
    
if ~isempty(tags) 
    nt = numel(tags);
    varargout = cell(1, nt);
    
    for i = 1 : numel(tags)
        
        t = tags(i);
        if t == 'C'
            varargout{i} = C;
        elseif t == 'I'
            varargout{i} = I;
        elseif t == 'G'
            varargout{i} = G;
        else
            error('uniqvalues:invalidarg', ...
                'Invalid output tag symbol.');                        
        end
    end    
end



