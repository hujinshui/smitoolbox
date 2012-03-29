function h = impline(a, b, c, rgn, varargin)
%IMPLINE Draw a line given by implicit equation
%
%   impline(a, b, c, [], ...);
%   impline(a, b, c, rgn, ...);
%
%       Draws a 2D line within a range. 
%       The line is given by the following implicit equation.
%
%           a x + b y + c = 0
%
%       In the input, rgn is the view range, in form of 
%       [left right top bottom].
%
%       rgn can also be input as [], then the range will be obtained 
%       from the current axis.
%
%       One can specify additional parameters, which will be forwarded
%       to the line function that this function actually invokes to 
%       do the drawing.
%

% Created by Dahua Lin, on Dec 31, 2011
%

%% verify inputs

if ~(isscalar(a) && isscalar(b) && isscalar(c))
    error('impline:invalidarg', 'a, b, and c should be all scalars.');
end

if isempty(rgn)
    rgn = [get(gca, 'XLim'), get(gca, 'YLim')];
else
    if ~(isnumeric(rgn) && isequal(size(rgn), [1 4]))
        error('impline:invalidarg', 'rgn should be a 1 x 4 numeric vector.');
    end
end

%% main

if abs(a) < abs(b)   
    x0 = rgn(1);
    x1 = rgn(2);
    y0 = - (a * x0 + c) / b;
    y1 = - (a * x1 + c) / b;    
else
    y0 = rgn(3);
    y1 = rgn(4);
    x0 = - (b * y0 + c) / a;
    x1 = - (b * y1 + c) / a;    
end

if nargout == 0
    line([x0 x1], [y0 y1], varargin{:});
else
    h = line([x0 x1], [y0 y1], varargin{:});
end



