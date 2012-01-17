function [theta, theta0] = mlogireg(X, y, w, rc, initsol, varargin)
%MLOGIREG Multi-class logistric regression
%
%   [theta, theta0] = MLOGIREG(X, y);
%   [theta, theta0] = MLOGIREG(X, y, w);
%   [theta, theta0] = MLOGIREG(X, y, w, rc, ...);
%   [theta, theta0] = MLOGIREG(X, y, w, rc, initsol, ...);
%
%       Performs multi-class logistic regression to find the discriminant
%       functions for classification.
%
%       Input arguments:
%       - X:        The feature matrix. Size: d x n, where d is the
%                   feature dimension, and n is the number of samples
%                   in X.
%
%       - y:        y can be in either of the following forms.
%                   - a label vector of size 1 x n. Here, y(i) can take
%                     value in {1, ..., K}. The function would set K
%                     to be max(y).
%                   - an assignment matrix of size K x n. Here, y(:,i)
%                     corresponds to X(:,i). Note that each column of y
%                     must sums to one.
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
%       - theta:    The coefficient matrix. Size: d x K
%       - theta0:   The vector of offsets. Size: K x 1.
%
%       The following expression yields a vector comprised of K scores,
%       with respect to the K classes.
%
%           theta' * x + theta0
%

% Created by Dahua Lin, on Jan 17, 2011
%

%% preprocess inputs

if nargin < 3
    w = [];
end
if nargin < 4
    rc = 0;
end
d = size(X, 1);

if size(y, 1) == 1
    K = max(y);
else
    K = size(y, 1);
end

if nargin < 5 || isempty(initsol)
    initsol = zeros((d+1) * K, 1);
else
    if ~(iscell(initsol) && numel(initsol) == 2)
        error('mlogireg:invalidarg', 'initisol is invalid for logireg.');
    end
    t = initsol{1};
    t0 = initsol{2};
    if ~(isfloat(t) && isequal(size(t), [d K]) && ...
            isfloat(t0) && isequal(size(t0), [K 1]))
        error('mlogireg:invalidarg', 'initisol is invalid for logireg.');
    end
    initsol = [t; t0.'];
    initsol = initsol(:);
end
           
opts = smi_optimset('bfgsfmin', varargin{:});

%% main

f = mlogireg_objfun(X, y, w, rc);
sol = bfgsfmin(f, initsol, opts);
sol = reshape(sol, d+1, K);

theta = sol(1:d, :);
theta0 = sol(d+1, :).';


