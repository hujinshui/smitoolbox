function L = wishartd_logpdf(Sigma, df, V, c, op)
%WISHARTD_LOGPDF Evaluates the log pdf of Wishart distribution
%
%   L = wishartd_logpdf(Sigma, df, V);
%   L = wishartd_logpdf(Sigma, df, V, 0);
%   L = wishartd_logpdf(Sigma, df, V, c);
%
%       Evaluates the log pdf at the matrices given in V, w.r.t. the
%       Wishart distribution.
%
%       Input arguments:
%       - Sigma:        the scale matrix (pdmat struct), Sigma.n == 1
%       - df:           the degree of freedom
%       - V:            the sample matrices (pdmat struct)
%       - c:            the constant term in log-pdf evaluation. If
%                       not supplied, the function will evaluate it.
%                       One can set c to zero to ignore this term.
%
%       Output arguments:
%       - L:            the result values in form of a row vector of
%                       size 1 x V.n.
%
%   L = wishartd_logpdf(InvSigma, df, V, [], 'inv');
%   L = wishartd_logpdf(InvSigma, df, V, 0, 'inv');
%   L = wishartd_logpdf(InvSigma, df, V, c, 'inv');
%
%       Here, InvSigma is the inverse scale matrix.
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify inputs

if ~(is_pdmat(Sigma) && Sigma.n == 1)
    error('wishartd_logpdf:invalidarg', ...
        'Sigma should be a pdmat struct with Sigma.d == 1.');
end
d = Sigma.d;

if ~(isfloat(df) && isscalar(df) && df > d - 1)
    error('wishartd_logpdf:invalidarg', ...
        'df should be a numeric scalar with df > d - 1.');
end

if ~(is_pdmat(V) && V.d == d)
    error('wishartd_logpdf:invalidarg', ...
        'V should be a pdmat struct with V.d == d.');
end

if nargin < 4 || isempty(c)
    calc_c = 1;
else
    if ~(isfloat(c) && isscalar(c) && isreal(c))
        error('wishartd_logpdf:invalidarg', 'c should be a real scalar.');
    end
    calc_c = 0;
end

if nargin < 5
    is_inv = 0;
else
    if ~(ischar(op) && strcmpi(op, 'inv'))
        error('wishartd_logpdf:invalidarg', ...
            'The 5th argument can only be ''inv''.');
    end
    is_inv = 1;
end

%% main

% inverse Sigma

if is_inv
    J = Sigma;
else
    J = pdmat_inv(Sigma);
end

% compute

L = (0.5) * ((df - d - 1) * pdmat_lndet(V) - pdmat_dot(J, V));

if calc_c
    c = wishartd_const(J, df);
end

if c ~= 0    
    L = L + c;    
end


