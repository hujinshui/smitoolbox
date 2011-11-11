function v = gentropy(C, op)
% Compute the entropy of a Gaussian distribution
%
%   v = gentropy(C);
%
%       computes the entropy of Gaussian distribution(s) based on
%       its covariance matrix.
%
%       C should be a pdmat struct. v will be a scalar (if C.n == 1),
%       or a 1 x C.n row vector.
%
%   v = gentropy(J, 'inv');
%
%       computes the entropy based on the information matrix.
%

% Created by Dahua Lin, on Sep 1, 2011
%

%% verify input

if ~is_pdmat(C)
    error('gentropy:invalidarg', 'The first arg should be a pdmat struct.');
end

if nargin < 2
    is_cov = 1;
else
    if ~(ischar(op) && strcmpi(op, 'inv'))
        error('gentropy:invalidarg', ...
            'The second arg to gentropy can only be ''inv''.');
    end
    is_cov = 0;
end

%% main

log2pip1 = 2.837877066409345483560659472811;

if is_cov            
    v = (C.d * log2pip1 + pdmat_lndet(C)) * 0.5;    
else
    v = (C.d * log2pip1 - pdmat_lndet(C)) * 0.5;    
end

