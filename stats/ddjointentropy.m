function v = ddjointentropy(cp, p0)
%DDJOINTENTROPY Computes the joint entropy for discrete distribution
%
%   v = ddjointentropy(cp, p0);
%       computes the joint entropy for distribution, based on the
%       conditional distribution p(y|x) and the prior distribution p0(x).
%
%       The joint entropy is defined as
%           
%           H(x, y) = - sum_{x,y} p(x,y) log p(x, y)
%
%       Suppose, there are m distinct values in x-space, and n distinct
%       values in y-space, then cp should be a m x n matrix, with cp(i, j)
%       equaling p(y = j | x = i), and thus all elements in sum(cp, 2)
%       should equal 1. And, p0 should be an m x 1 vector, with sum(p0)
%       being 1.
%
%       The output is the joint entropy of the distribution characterized
%       by the prior and conditional distributions.
%
%   v = ddjointentropy(cp);
%       computes the joint entropy based on the conditional probability
%       given in cp, and equal prior.
%

%   History
%       - Created by Dahua Lin, on Jun 5, 2008
%

%% parse and verify input arguments

assert(isfloat(cp) && ndims(cp) == 2, 'ddjointentropy:invalidarg', ...
    'cp should be a numeric matrix.');

m = size(cp, 1);

if nargin < 2 || isempty(p0)    
    equal_prior = true;
    
else
    assert(isfloat(p0) && ndims(p0) == 2 && size(p0,1) == m && size(p0,2) == 1, ...
        'ddjointentropy:invalidarg', ...
        'p0 should be a m x 1 numeric vector.');
    
    equal_prior = false;    
end

%% main

% compute joint distribution
if equal_prior    
    jp = cp * (1 / m);
    
else
    jp = bsxfun(@times, cp, p0);
end

% compute joint entropy
v = ddentropy(jp(:));