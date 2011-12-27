function Gpos = gaussd_conjupdate(Gpri, S)
% Computes the posterior Gaussian distribution via conjuate update
%
%   Gpos = gaussd_conjupdate(Gpri, S);
%
%       Estimates the posterior distribution via conjugate updates on
%       a given priro.
%       
%       Input arguments:
%       - G:    the Gaussian prior (G.ty == 'c' && G.n == 1)
%       - S:    the gaussd object that captures the updates
%
%       Output arguments:
%       - Gpos:     the posterior Gaussian distribution
%       

% Created by Dahua Lin, on Dec 14, 2011
%

%% verify input arguments

if ~(is_gaussd(Gpri) && Gpri.ty == 'c' && Gpri.n == 1)
    error('gaussd_conjupdate:invalidarg', ...
        'G should be a gaussd struct with G.ty == ''c'' and G.n == 1.');
end

d = Gpri.d;

if ~(is_gaussd(S) && S.d == d && S.ty == 'c')
    error('gaussd_conjupdate:invalidarg', ...
        'S should be a gaussd struct with G.ty == ''c'' and G.n == 1.');
end


%% main

% posterior h

dh = S.h;
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

J = pdmat_plus(Gpri.J, S.J);

% output

Gpos.tag = 'gaussd';
Gpos.ty = 'c';
Gpos.n = size(h, 2);
Gpos.d = d;
Gpos.h = h;
Gpos.J = J;


