function v = Lpnorm(X, p, d)
%Compute Lp-norm values of vectors
%
%   v = Lpnorm(X, p);
%   v = Lpnorm(X, p, d);
%       compute the Lp-norm values of the vectors contained in X along
%       dimension d. If d is omitted, by default, it is set to 1.
%
%   v = Lpnorm(X, p, []);
%       compute the Lp-norm of the entire array X.
%
%   Remarks:
%       - For this function, p should be a scalar that has p >= 1.
%

% Created by Dahua Lin, on Aug 1, 2010
%

%% verify input

if nargin < 3
    d = 1;
end

if p < 1
    error('Lpnorm:invalidarg', 'p must have p >= 1.');
end

%% main

if p == 1
    v = L1norm(X, d);
    
elseif p == 2
    v = L2norm(X, d);
    
elseif isinf(p)
    v = Linfnorm(X, d);
    
else
    if ~isempty(d)
        v = sum(abs(X).^p, d);
    else
        v = sum(abs(X(:)).^p);
    end

    v = v.^(1/p);
end

