function v = gaussd_entropy(C, op)
% Compute the entropy of a Gaussian distribution
%
%   v = gaussd_entropy(C);
%   v = gaussd_entropy(C, 'm');
%
%       computes the entropy of Gaussian distribution(s) based on
%       its covariance matrix.
%
%       C should be a pdmat struct. v will be a scalar (if C.n == 1),
%       or a 1 x C.n row vector.
%
%   v = gaussd_entropy(J, 'i');
%   v = gaussd_entropy(J, 'c');
%
%       computes the entropy based on the information matrix.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 1, 2011
%       - Modified by Dahua Lin, on Dec 5, 2011
%

%% verify input

if ~is_pdmat(C)
    error('gaussd_entropy:invalidarg', 'The first arg should be a pdmat struct.');
end

if nargin < 2
    is_cov = 1;
else
    if ~(ischar(op) && isscalar(op))
        error('gaussd_entropy:invalidarg', ...
            'The second arg to gauss_entropy should be a character.');
    end
    if op == 'm'
        is_cov = 1;
    elseif op == 'i' || op == 'c'
        is_cov = 0;
    else
        error('gaussd_entropy:invalidarg', 'The 2nd arg to gauss_entropy is invalid.');
    end
end

%% main

log2pip1 = 2.837877066409345483560659472811;

if is_cov            
    v = (C.d * log2pip1 + pdmat_lndet(C)) * 0.5;    
else
    v = (C.d * log2pip1 - pdmat_lndet(C)) * 0.5;    
end

