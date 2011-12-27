function c = wishartd_const(J, df)
% Compute the constant term for Wishart distribution log-pdf
%
%   c = wishartd_const(J, df);
%
%       Evaluates the constant term for Wishart distribution log-pdf.
%
%       Inputs:
%       - J:        inverse scale matrix in form of pdmat struct.
%       - df:       the degree of freedom.
%

% Created by Dahua Lin, on Dec 26, 2011
%

%% verify input

if ~(is_pdmat(J) && J.n == 1)
    error('wishartd_const:invalidarg', ...
        'J should be a pdmat struct with J.n == 1.');
end
d = J.d;

if ~(isfloat(df) && isreal(df) && isscalar(df) && df > d - 1)
    error('wishartd_const:invalidarg', ...
        'df should be a real scalar greater than d - 1.');
end

%% main

c = (df / 2) * pdmat_lndet(J) - (df * d / 2) * log(2) - mvgammaln(d, df / 2);

