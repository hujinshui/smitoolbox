function [h, J, u] = gaussgm_pos(gpri, X, Jx, A, Z)
% Compute the posterior distribution of a Gaussian generative model
%
%   The Gaussian generative model is formulated as follows
%
%       u ~ Gauss(mu0, Cu);
%
%       x ~ Gauss(u, Cx);       -- Formulation (1)
%    or x ~ Gauss(A * u, Cx);   -- Formulation (2)
%
%
%   [h, J] = gaussgm_pos(gpri, X, Jx);      -- for Formulation (1)
%   [h, J] = gaussgm_pos(gpri, X, Jx, A);   -- for Formulation (2)
%
%       returns the potential vector and the information matrix of 
%       the posterior distribution of u, given a collection of 
%       observed samples, as columns of X.
%
%       Input arguments:
%       - gpri:     the prior distribution, which should be a gaussd
%                   object. 
%                   One can set gpri to [], to use an uninformative prior.
%
%       - X:        the observed sample matrix, whose size is d x n.
%                   Here, d is the space dimension, and n is the number
%                   of samples. Each column of X is a sample.
%
%       - Jx:       the inverse covariance (pdmat struct) the conditional
%                   distribution. Jx = inv(Cx).
%
%       - A:        the transform matrix (for Formulation (2)), or 
%                   a scalar.
%
%       Output arguments:
%       - h:        the posterior potential vector
%
%       - J:        the posterior information matrix 
%
%
%   [h, J] = gaussgm_pos(gpri, X, Jx, [], w);
%   [h, J] = gaussgm_pos(gpri, X, Jx, A, w);
%
%       computes the posterior parameters with weighted samples. 
%       w is either a 1 x n row vector or a K x n matrix that gives
%       K different sets of weights.
%
%       Correspondingly, h is either a d x 1 vector of a d x K matrix,
%       J.n = 1 or K. 
%
%   [h, J] = gaussgm_pos(gpri, X, Jx, [], g);
%   [h, J] = gaussgm_pos(gpri, X, Jx, A, g);
%
%       computes the posterior parameters with grouped samples. 
%       In particular, the k-th column of h and the J-th matrix in J
%       corresponds to the posterior distribution conditioned on the
%       samples selected by g{k}.
%
%       Here, g is a cell array of index vectors.
%
%   [h, J, u] = gaussgm_pos( ... );
%
%       additionally returns the MAP estimates of the parameter u.
%

% Created by Dahua Lin, on Aug 26, 2011
%
 
%% verify input arguments

if ~isempty(gpri)
    if ~(isa(gpri, 'gaussd') && gpri.num == 1 && gpri.has_ip)
        error('gaussgm_pos:invalidarg', ...
            'gpri should be an instance of class gaussd with gpri.num == 1.');
    end
    has_pri = true;
else
    has_pri = false;    
end

if ~(isfloat(X) && ndims(X) == 2 && isreal(X))
    error('gaussgm_pos:invalidarg', 'X should be a real matrix.');
end
[xdim, n] = size(X);

if ~(is_pdmat(Jx) && Jx.n == 1)
    error('gaussgm_pos:invalidarg', ...
        'Jx should be a pdmat struct with single matrix.');
end

if Jx.d ~= xdim
    error('gaussgm_pos:invalidarg', ...
        'The dimension of X is not consistent with Jx.d.');
end


if nargin < 4
    A = [];
end

if isempty(A) || isscalar(A)   
    if gpri.dim ~= xdim
        error('gaussgm_pos:invalidarg', ...
            'gpri.dim is inconsistent with the sample dimension.');
    end
else
    if ~( isfloat(A) && isreal(A) && ndims(A) == 2 && size(A,1) == xdim ) 
        error('gaussgm_pos:invalidarg', ...
            'A should be either a scalar or a real matrix with size(A,1) == size(X,1).');
    end
    
    if has_pri
        if size(A, 2) ~= gpri.dim
            error('gaussgm_pos:invalidarg', ...
                'The size of A is inconsistent with gpri.dim.');
        end
    end
end

if nargin < 5
    zty = 0;   % no grouping and weighting
    K = 1;
else
    if isnumeric(Z)
        zty = 1;    % weighting
        w = Z;
        if ~(isfloat(w) && ndims(w) == 2 && isreal(w) && size(w,2) == n)
            error('gaussgm_pos:invalidarg', ...
                'w should be a real matrix with size(w,2) == n.');
        end
        K = size(w, 1);
        
    elseif iscell(Z)
        zty = 2;    % grouping
        g = Z;
        K = numel(g);
        
    end
end

%% main

% formulas to implement:
%
% h_k <- A' * Jx * (sum_i w_{ki} x_i)
% J_k <- (sum_i w_i) A' * Jx * A
%

% compute sample sum and weight sum

if zty == 0
    sx = sum(X, 2);
    sw = n;
    
elseif zty == 1  % weighting
    sx = X * w';
    sw = sum(w, 2).';
    
else  % grouping
    sx = zeros(xdim, K);
    sw = zeros(1, K);
    for k = 1 : K
        Xk = X(:, g{k});
        sx(:, k) = sum(Xk, 2);
        sw(k) = size(Xk, 2);
    end
    
end

% compute h and J (based on samples)
 
h = pdmat_mvmul(Jx, sx);
if ~isempty(A)
    h = A' * h;
end

if isempty(A) 
    J = pdmat_scale(Jx, sw);
    
elseif isscalar(A)
    J = pdmat_scale(Jx, sw * (A^2));
    
else  % A is a transform matrix
    J = pdmat_quad(Jx, A);
    J = pdmat_scale(pdmat(J), sw);    
end

% incorporate prior

if ~isempty(gpri)
    
    if K == 1
        h = gpri.h + h;
    else
        h = bsxfun(@plus, gpri.h, h);
    end
       
    J = pdmat_plus(gpri.J, J);    
end
 
% MAP estimation of u

if nargout >= 3    
    u = pdmat_lsolve(J, h);
end


