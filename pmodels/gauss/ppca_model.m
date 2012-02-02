function M = ppca_model(B, s, se, mu)
% Constructs a Probabilistic PCA model struct
%
%   M = ppca_model(B, s, se);
%   M = ppca_model(B, s, se, mu);
%
%       constructs a PPCA model strict.
%
%       Input arguments:
%
%       - B:        the basis matrix. The size of U is d x q, where
%                   d is the observed space dimension, and q is the
%                   latent space dimension.
%
%                   Note that B should be an orthonormal matrix, meaning
%                   that B' * B is identity.
%
%       - s:        the vector of principal scales, i.e. the scale along 
%                   principal directions. The length of s is q.
%
%                   Note that the factor matrix W has W = B * diag(s).
%
%       - se:       the noise standard deviation.
%   
%       - mu:       the mean vector. It can either a d x 1 vector or
%                   a zero scalar. (If omitted, mu is set to zero).

%   History
%   -------
%       - Created by Dahua Lin, on Nov 20, 2010
%       - Modified by Dahua Lin, on Nov 3, 2011
%       - Modified by Dahua Lin, on Dec 27, 2011
%

%% verify inputs

if ~(isfloat(B) && isreal(B) && ndims(B) == 2)
    error('ppca_model:invalidarg', 'U should be a real matrix.');
end
[d, q] = size(B);

if q >= d
    error('ppca_model:invalidarg', 'q should be less than d.');
end

if ~(isfloat(s) && isreal(s) && isvector(s) && numel(s) == q)
    error('ppca_model:invalidarg', 's should be a vector of length q.');
end
if size(s, 1) > 1
    s = s.';    % turn into row vector
end

if ~(isfloat(se) && isreal(se) && se > 0)
    error('ppca_model:invalidarg', 'se should be a positive real scalar.');
end

if nargin < 4
    mu = 0;
else
    if ~(isequal(mu, 0) || ...
            (isfloat(mu) && isreal(mu) && isequal(size(mu), [d 1])))
        error('ppca_model:invalidarg', ...
            'mu should be either zero or a d x 1 real vector.');
    end
end


%% main

% basic fields

M.tag = 'ppca';
M.d = d;
M.q = q;
M.mu = mu;
M.B = B;
M.s = s;
M.se = se;

% auxiliary fields

evs = zeros(1, d);
evs(1:q) = s.^2 + se^2;
evs(q+1:end) = se^2;
M.ldc = sum(log(evs));

