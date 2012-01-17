function [theta, theta0] = logireg(X, y, w, rc, initsol, varargin)
%LOGIREG Logistic regression
%
%   [theta, theta0] = LOGIREG(X, y);
%   [theta, theta0] = LOGIREG(X, y, w);
%   [theta, theta0] = LOGIREG(X, y, w, rc, ...);
%   [theta, theta0] = LOGIREG(X, y, w, rc, initsol, ...);
%
%       Performs logistic regression to find a decision boundary
%       between two classes. 
%
%       Input arguments:
%       - X:        The feature matrix. Size: d x n, where d is the
%                   feature dimension, and n is the number of samples
%                   in X.
%
%       - y:        The vector of class labels, whose values can be
%                   either 1 or -1. The length of y should be n.
%
%       - w:        The sample weights. It is either empty (indicating
%                   that all samples have the same weight), or a vector
%                   of length n that gives sample-specific weights.
%
%       - rc:       The regularization coefficient. 
%                   (If omitted, rc is set to zero).
%
%       - initsol:  The initial solution in form of {theta, theta0}.
%                   This can be omitted or input as an empty array, in
%                   which case, both theta and theta0 are initialized
%                   to be zeros.
%
%       One can further specify the following options in form of 
%       name/value pairs to control the optimization procedure. 
%      
%       - MaxIter:  the maximum number of iterations {100}
%       - TolFun:   the termination tolerance of objective value change {1e-8}
%       - TolX:     the termination tolerance of solution change {1e-7}
%       - Display:  the level of display {'none'}|'proc'|'iter'
%
%       Output arguments:
%       - theta:    The coefficient vector. Size: d x 1
%       - theta0:   The offset scalar.
%
%       The linear predictor is then given by 
%
%           theta' * x + theta0
%
%       One can determine which class x belongs to based on the sign
%       of this predictor.
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% pre-process inputs

if nargin < 3
    w = [];
end
if nargin < 4
    rc = 0;
end
d = size(X, 1);
if nargin < 5 || isempty(initsol)
    initsol = zeros(d+1, 1);
else
    if ~(iscell(initsol) && numel(initsol) == 2)
        error('logireg:invalidarg', 'initisol is invalid for logireg.');
    end
    t = initsol{1};
    t0 = initsol{2};
    if ~(isfloat(t) && isequal(size(t), [d 1]) && ...
            isfloat(t0) && isscalar(t0))
        error('logireg:invalidarg', 'initisol is invalid for logireg.');
    end
    initsol = [t; t0];
end
           
opts = smi_optimset('bfgsfmin', varargin{:});

%% main

f = logireg_objfun(X, y, w, rc);
sol = bfgsfmin(f, initsol, opts);
theta = sol(1:d);
theta0 = sol(d+1);


