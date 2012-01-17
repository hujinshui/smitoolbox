function [theta, theta0] = glinreg(X, Y, w, rho, rc, initsol, varargin)
%GLINREG Logistic regression
%
%   [theta, theta0] = GLINREG(X, Y, [], rho);
%   [theta, theta0] = GLINREG(X, Y, w, rho);
%   [theta, theta0] = GLINREG(X, Y, w, rho, rc, ...);
%   [theta, theta0] = GLINREG(X, Y, w, rho, rc, initsol, ...);
%
%       Performs generalized linear regression.
%
%       Input arguments:
%       - X:        The feature matrix. Size: d x n, where d is the
%                   feature dimension, and n is the number of samples
%                   in X.
%
%       - Y:        The matrix of output vectors. Size: q x n, where q
%                   is the dimension of the output space.
%
%       - w:        The sample weights. It is either empty (indicating
%                   that all samples have the same weight), or a vector
%                   of length n that gives sample-specific weights.
%
%       - rho:      The loss function handle, which takes difference
%                   vectors as input, and outputs the loss values.
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
%       - theta:    The coefficient vector. Size: d x q
%       - theta0:   The offset scalar. Size: q x 1.
%
%       The predicted output for a given input x is
%
%           theta' * x + theta0
%

% Created by Dahua Lin, on Jan 17, 2012
%

%% preprocess arguments

if nargin < 3
    w = [];
end
if nargin < 5
    rc = 0;
end
d = size(X, 1);
q = size(Y, 1);

if nargin < 6 || isempty(initsol)
    initsol = zeros((d+1) * q, 1);
else
    if ~(iscell(initsol) && numel(initsol) == 2)
        error('glinreg:invalidarg', 'initisol is invalid for logireg.');
    end
    t = initsol{1};
    t0 = initsol{2};
    if ~(isfloat(t) && isequal(size(t), [d q]) && ...
            isfloat(t0) && isequal(size(t0), [q 1]))
        error('glinreg:invalidarg', 'initisol is invalid for logireg.');
    end
    initsol = [t; t0.'];
    initsol = initsol(:);
end
           
opts = smi_optimset('bfgsfmin', varargin{:});

%% main

f = glinreg_objfun(X, Y, w, rho, rc);
sol = bfgsfmin(f, initsol, opts);
theta = sol(1:d, :);
theta0 = sol(d+1, :).';


