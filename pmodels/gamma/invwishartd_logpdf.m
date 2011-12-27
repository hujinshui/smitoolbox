function L = invwishartd_logpdf(Phi, df, V, c)
%INVWISHARTD_LOGPDF Evaluates the log pdf of Inverse Wishart distribution
%
%   L = invwishartd_logpdf(Phi, df, V);
%   L = invwishartd_logpdf(Phi, df, V, 0);
%   L = invwishartd_logpdf(Phi, df, V, c);
%
%       Evaluates the log pdf at the matrices given in V, w.r.t. the
%       inverse Wishart distribution.
%
%       Input arguments:
%       - Phi:          the inverse scale matrix (pdmat struct), Phi.n == 1
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

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify inputs

if ~(is_pdmat(Phi) && Phi.n == 1)
    error('invwishartd_logpdf:invalidarg', ...
        'Phi should be a pdmat struct with Sigma.d == 1.');
end
d = Sigma.d;

if ~(isfloat(df) && isscalar(df) && df > d - 1)
    error('invwishartd_logpdf:invalidarg', ...
        'df should be a numeric scalar with df > d - 1.');
end

if ~(is_pdmat(V) && V.d == d)
    error('invwishartd_logpdf:invalidarg', ...
        'V should be a pdmat struct with V.d == d.');
end

if nargin < 4 || isempty(c)
    calc_c = 1;
else
    if ~(isfloat(c) && isscalar(c) && isreal(c))
        error('invwishartd_logpdf:invalidarg', 'c should be a real scalar.');
    end
    calc_c = 0;
end

%% main

% compute

IV = pdmat_inv(V);
L = (-0.5) * ((df + d + 1) * pdmat_lndet(V) - pdmat_dot(Phi, IV));

if calc_c
    c = wishartd_const(Phi, df);
end

if c ~= 0    
    L = L + c;    
end


