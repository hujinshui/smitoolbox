function Pi = fmm_estPi(K, Z, alpha, op)
% Estimates the component Prior given sample labels
%
%   Pi = fmm_estPi(K, Z);
%   Pi = fmm_estPi(K, Z, alpha);
%       
%       Estimates the prior probabilities of the mixture components (Pi).
%       Given:
%
%       - K:        the number of components
%       - Z:        the labeling, which is either a 1 x n label vector,
%                   or a K x n soft assignment matrix.
%       - alpha:    the concentration parameter of Dirichlet prior
%                   (by default, it is 1)  
%
%       Note that alpha can be inf, in which case, Pi is always a 
%       uniform distribution.
%
%   Pi = fmm_estPi( ..., 'sample');
%
%       Samples Pi from the posterior Dirichlet distribution.
%              

% Created by Dahua Lin, on Sep 28, 2011
%

%% verify input arguments

if ~(isnumeric(K) && isscalar(K) && K == fix(K) && K >= 1)
    error('fmm_estPi:invalidarg', 'K should be a positive integer scalar.');
end

if ~(isnumeric(Z) && isreal(Z) && ndims(Z) == 2)
    error('fmm_estPi:invalidarg', 'Z should be a real vector or matrix.');
end

if nargin < 3
    alpha = 1;
else
    if ~(isfloat(alpha) && isreal(alpha) && isscalar(alpha) && alpha >= 1)
        error('fmm_estPi:invalidarg', ...
            'alpha should be a real value in [1, inf].');
    end
end

if nargin < 4
    samp = false;
else
    if ~(ischar(op) && strcmpi(op, 'sample'))
        error('fmm_estPi:invalidarg', ...
            'The 4th argument to fmm_estPi can only be ''sample''.');
    end
    samp = true;
end
   

%% main

if K == 1
    Pi = 1;
    return;
end

if isinf(alpha)
    Pi = constmat(K, 1, 1 / double(K));
    return;
end

if size(Z, 1) == 1  % label vector
    tw = intcount(K, Z).';    
elseif size(Z, 1) == K  % soft assignment matrix
    tw = sum(Z, 2);    
else
    error('fmm_estPi:invalidarg', 'The size of Z is invalid.');
end
    
if ~samp
    if alpha ~= 1
        tw = tw + (alpha - 1);
    end
    Pi = tw ./ sum(tw);    
else
    Pi = dirichlet_sample(K, tw + alpha, 1);
end


