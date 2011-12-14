function varargout = gaussd_conjupdate(Gpri, dh, dJ)
% Computes the posterior Gaussian distribution via conjuate update
%
%   Gpos         = gaussd_conjupdate(Gpri, dh, dJ);
%   [hpos, Jpos] = gaussd_conjupdate(Gpri, dh, dJ);
%
%       Estimates the posterior distribution via conjugate updates on
%       a given priro.
%       
%       Input arguments:
%       - G:    the Gaussian prior (G.ty == 'c' && G.n == 1)
%       - dh:   the update to the potential vector
%       - dJ:   the update to the precision matrix.
%
%       Output arguments:
%       - Gpos:     the posterior Gaussian distribution
%       - hpos:     the potential vector of the posterior Gaussian(s)
%       - Jpos:     the precision matrix of the posterior Gaussian(s)
%       

% Created by Dahua Lin, on Dec 14, 2011
%

%% verify input arguments

if ~(is_gaussd(Gpri) && Gpri.ty == 'c' && Gpri.n == 1)
    error('gaussd_conjupdate:invalidarg', ...
        'G should be a gaussd struct with G.ty == ''c'' and G.n == 1.');
end

d = Gpri.d;
if ~(isfloat(dh) && isreal(dh) && ndims(dh) == 2 && size(dh, 1) == d)
    error('gaussd_conjupdate:invalidarg', ...
        'dh should be a real matrix with G.d rows.');
end

if ~(is_pdmat(dJ) && dJ.d == d)
    error('gaussd_conjupdate:invalidarg', ...
        'dJ should be a pdmat struct with J.d == d.');
end

n = size(dh, 2);
if ~(dJ.n == 1 || dJ.n == n)
    error('gaussd_conjupdate:invalidarg', ...
        'dJ.n and the number of columns in dh are not consistent.');
end

%% main

% posterior h

if isequal(Gpri.h, 0)
    h = dh;
else
    if n == 1
        h = Gpri.h + dh;
    else
        h = bsxfun(@plus, Gpri.h, dh);
    end
end

% posterior J

J = pdmat_plus(Gpri.J, dJ);

% output

if nargout <= 1
    Gpos = gaussd('c', h, J);
    varargout = {Gpos};
else
    varargout = {h, J};
end

