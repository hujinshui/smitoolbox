function fc = combfun(funs, w)
% Combining multiple objective functions into one
%
%   fc = combfun({f1, f2, ...});
%       The first argument to this function is a cell array containing
%       objective function handles.
%
%       This statement returns a combined function handle fc, which equals
%       f1 + f2 + ...
%
%       Note that if each of the function fi yields n output arguments,
%       then fc yields n output arguments correspondingly.
%
%   fc = combfun({f1, f2, ...}, [c1, c2, ...]);
%       makes a weighted combination of the supplied functions using 
%       the weights in the second input argument, which is a vector.
%
%       The resultant function is
%
%           f = c1 * f1 + c2 * f2 + ...
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 24, 2011
%

%% main

if ~(iscell(funs) && ~isempty(funs))
    error('combfun:invalidarg', ...
        'The 1st arg to combfun should be a cell array of function handles.');
end

n = numel(funs);
for i = 1 : n
    if ~isa(funs{i}, 'function_handle')
        error('combfun:invalidarg', ...
            'The 1st arg to combfun should be a cell array of function handles.');
    end
end

if nargin < 2
    w = ones(1, n);
else
    if ~(isfloat(w) && isvector(w) && numel(w) == n)
        error('combfun:invalidarg', 'The weight vector is invalid.');
    end
end

%% main

fc = @(x) fcomb(x, funs, w);

%% The combined function

function varargout = fcomb(x, funs, w)

nout = max(nargout, 1);
n = numel(funs);

f1 = funs{1};
[vout{1:nout}] = f1(x);
if w(1) ~= 1
    for j = 1 : nout
        vout{j} = w(1) * vout{j};
    end
end
for i = 2 : n
    [cv{1:nout}] = funs{i}(x); %#ok<AGROW>
    if w(i) == 1
        for j = 1 : nout
            vout{j} = vout{j} + cv{j};
        end
    else
        for j = 1 : nout
            vout{j} = vout{j} + w(i) * cv{j};
        end
    end
end

varargout = vout;

     
    