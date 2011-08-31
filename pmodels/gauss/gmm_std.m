function prg = gmm_std(X, K, Cx, gpri, iwspri, method)
% Creates a standard Gaussian mixture model (GMM) program
%
%   The GMM program is an SMI program that implements the inference
%   over the following formulation.
%
%       u_k ~ Gauss(mu, Cu);        for k = 1, ..., K
%       x_i ~ Gauss(u_{z_i}, Cx);   for i = 1, ..., N
%
%   Cx is either given, or from an inverse Wishart distribution.
%
%   prg = gmm_std(X, K);
%   prg = gmm_std(X, K, Cx);
%   prg = gmm_std(X, K, [], gpri);
%   prg = gmm_std(X, K, Cx, gpri);
%
%   prg = gmm_std(X, K, [], gpri, iwspri); (not implemented yet)
%
%   prg = gmm_std( ......, method);
%       
%       Construsts an SMI program that implements the inference of
%       GMM parameters. 
%       
%       Input arguments:
%       - X:        the observed data [xdim x n]
%       - K:        the number of mixture components
%       - Cx:       the covariance for generating x.
%                   (If left empty or omitted, it is to be inferred)
%
%       - gpri:     the Gaussian prior for the mean
%       - iwspri:   the Inverse Wishart prior for the covariance
%
%       - method:   the inference method:
%                   - 'em': using expectation-maximization (default)
%                   -' gs': using gibbs sampling
%

% Created by Dahua Lin, on Aug 31, 2011
%

%% verify input arguments

if ~(isfloat(X) && ndims(X) == 2 && ~issparse(X) && isreal(X))
    error('gmm_std:invalidarg', 'X should be a real matrix.');
end
[d, N] = size(X);

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('gmm_std:invalidarg', 'K should be a non-negative integer scalar.');
end

if nargin < 3 || isempty(Cx)
    est_Cx = true;
    Cx = [];
else
    est_Cx = false;
    if ~(ispdmat(Cx) && Cx.d == d)
        error('gmm_std:invalidarg', ...
            'Cx should be a pdmat struct with Cx.dim == size(X,1).');
    end
end

if nargin < 4 || isempty(gpri)
    mu = [];
    Cu = [];
else
    if ~(isa(gpri, 'gaussd') && gpri.has_mp)
        error('gmm_std:invalidarg', ...
            'gpri should be a gaussd object with mean parameters.');
    end
    if ~(gpri.num == 1 && gpri.dim == d)
        error('gmm_std:invalidarg', 'gpri is invalid.');
    end
    mu = gpri.mu;
    Cu = gpri.C;
end
       
if nargin < 6
    method = 'em';
else
    if ~isvarname(method)
        error('gmm_std:invalidarg', 'method is not a valid name.');
    end
    method = lower(method);
    if ~(strcmp(method, 'em') || strcmp(method, 'gs'))
        error('gmm_std:invalidarg', 'Unknown method name %s', method);
    end
end

if est_Cx
    error('gmm_std:notimplement', ...
        'Estimation of Cx has not yet been implemented.');
end


%% main

% create model

gm = gaussgm(Cx, Cu);


    
    
    


