function clsname = fresult_class(a, b, c)
% determine the class of the floating point computation result 
%
%   clsname = fresult_class(a, b);
%       determines the result class of a * b
%
%   clsname = fresult_class(a, b, c);
%       determines the result class of a * b * c
%

% Created by Dahua Lin, on Aug 14, 2011
%

%% verify input

if ~isfloat(a)
    error('fresult_class:invalidarg', 'a should be of floating-point type.');
end

if ~isfloat(b)
    error('fresult_class:invalidarg', 'b should be of floating-point type.');
end

if nargin >= 3
    if ~isfloat(c)
        error('fresult_class:invalidarg', 'c should be of floating-point type.');
    end
end

%% main

if nargin < 3    
    if isa(a, 'double') && isa(b, 'double')
        clsname = 'double';
    else
        clsname = 'single';
    end
    
else    
    if isa(a, 'double') && isa(b, 'double') && isa(c, 'double')
        clsname = 'double';
    else
        clsname = 'single';
    end
end

    


