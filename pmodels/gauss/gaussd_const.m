function [ca, cb] = gaussd_const(G)
% Calculate two useful constants of Gaussian distributions
%
%   ca = gaussd_const(G);
%   [ca, cb] = gaussd_const(G);
%
%       This function calculates two useful constants for Gaussian
%       pdf evaluation.
%
%       ca = mu' * h = mu' * inv(C) * mu = h' * inv(J) * h
%       cb = -(1/2) * log((2 * pi)^d * |C|) = d/2 - entropy
%
%       If there is only one output arguments, then only ca is computed.
%
%       Outputs:
%       - ca:       In general, it will be a vector of size 1 x G.n
%                   If G.mu or G.h is a zero scalar, then ca is a 
%                   zero scalar.
%
%       - cb:       It is a scalar if n' is 1, or a vector of size 1 x n'.
%                   Here n' is G.C.n or G.J.n.
%   

% Created by Dahua Lin, on Dec 5, 2011
%

%% verify

if ~is_gaussd(G)
    error('gaussd_const:invalidarg', 'G must be a gaussd struct.');
end

%% main

ty = G.ty;

% calculate ca

if ty == 'm'
    ca = calc_ca(G.n, G.mu, G.C);
else
    ca = calc_ca(G.n, G.h, G.J);
end
        
% calculate cb

if nargout >= 2
    if ty == 'm'
        cb = G.d * log(2 * pi) + pdmat_lndet(G.C);
    else
        cb = G.d * log(2 * pi) - pdmat_lndet(G.J);
    end
    cb = (-0.5) * cb;
end


%% core functions

function ca = calc_ca(n, u, S)

if isequal(u, 0)
    ca = 0;
else
    if n == 1
        ca = u' * pdmat_lsolve(S, u);
    else
        h = pdmat_lsolve(S, u);
        ca = dot(u, h, 1);
    end
end

