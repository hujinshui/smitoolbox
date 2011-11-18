function f = robustloss(name)
% Get a robust loss function
%
%   f = robustloss(name);
%       Get the specified robust loss function, which converts a residue
%       value to a loss value. 
%
%       It returns a function handle f, which supports the following
%       usages:
%           
%           v = f(r);
%           [v, dv1] = f(r);
%           [v, dv1, dv2] = f(r);
%           w = f(r, 'w');
%
%       f takes into the residue r as input, and outputs the function
%       value v, 1st-order derivative dv1, and the 2nd-order derivative
%       dv2.
%
%       When the 2nd argument is 'w', it returns the weight value,
%       which equals dv1 ./ v.
%
%
%   The following is a list of the names of available loss functions:
%
%       - 'bisquare':   f(r) = (1 - (1 - r^2)^3) / 6
%                       w(r) = (1 - r^2)^2 (|r| <= 1) or 0 (|r| > 1)
%
%       - 'huber':      f(r) = r^2 / 2 (|r| <= 1) or |r| - 1/2 (|r| > 1)
%                       w(r) = 1 / max(1, |r|)
%
%       - 'fair':       f(r) = |r| - log(1 + |r|)
%                       w(r) = 1 / (1 + |r|)
%
%       - 'logistic':   f(r) = log(cosh(r))
%                       w(r) = tanh(r) / r
%
%       - 'talwar':     f(r) = min(r, 1)^2 / 2
%                       w(r) = 1 (r <= 1) or 0 (r > 1)
%

% Created by Dahua Lin, on Jan 5, 2011
%

%% verify input

if ~ischar(name)
    error('robustloss:invalidarg', 'The 1st arg should be a char string.');
end


%% main

switch lower(name)
    case 'bisquare'
        f = @bisquare;
        
    case 'huber'
        f = @huber;
        
    case 'fair'
        f = @fair;
        
    case 'logistic'
        f = @logistic;
        
    case 'talwar'
        f = @talwar;
        
    otherwise
        error('robustloss:invalidarg', ...
            'Unknown robust loss function name %s', name);
end


%% robust loss functions

% bisquare

function [v, dv1, dv2] = bisquare(r, op)

r2 = r.^2;
i = abs(r) <= 1;

if nargin == 1
    v = (1/6) * (1 - ((1 - r2).^3) .* i);

    if nargout >= 2
        dv1 = r .* (1 - r2).^2 .* i;
    end
    if nargout >= 3
        dv2 = (1 - 6 * r2 + 5 * r2.^2) .* i;
    end
    
elseif isequal(op, 'w')
    v = ((1 - r2).^2) .* i;
    
end


% huber

function [v, dv1, dv2] = huber(r, op)

ar = abs(r);

if nargin == 1
    b = min(ar, 1);
    v = b .* (ar - 0.5 * b);
    if nargout >= 2
        dv1 = sign(r) .* b;
    end
    if nargout >= 3
        dv2 = ar <= 1;
    end
    
elseif isequal(op, 'w')
    v = 1 ./ max(1, ar);
    
end



% fair

function [v, dv1, dv2] = fair(r, op)

ar = abs(r);

if nargin == 1
    v = ar - log(1 + ar);
    if nargout >= 2
        dv1 = r ./ (1 + ar);
    end
    if nargout >= 3
        dv2 = 1 ./ (1 + ar).^2;
    end
    
elseif isequal(op, 'w')
    v = 1 ./ (1 + ar);
    
end


% logistic

function [v, dv1, dv2] = logistic(r, op)

if nargin == 1
    cr = cosh(r);
    v = log(cr);    
    if nargout >= 2
        dv1 = tanh(r);
    end
    if nargout >= 3
        dv2 = cr .^ (-2);
    end
    
elseif isequal(op, 'w')
    v = tanh(r) ./ r;
    
end
    

% talwar

function [v, dv1, dv2] = talwar(r, op)

ar = abs(r);

if nargin == 1
    i = ar <= 1;
    v = 0.5 * (min(ar, 1).^2);
    if nargout >= 2
        dv1 = r .* i;
    end
    if nargout >= 3
        dv2 = double(i);
    end
    
elseif isequal(op, 'w')
    v = double(ar <= 1);
    
end





