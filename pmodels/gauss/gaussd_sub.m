function Gs = gaussd_sub(G, idx)
% Get a subset of Gaussian models
%
%   Gs = gaussd_sub(G, idx);
%       Get a subset of Gaussian models (specified by the index or 
%       index vector idx).
%

% Created by Dahua Lin, on Dec 20, 2011
%

%% main

if G.ty == 'm'
    mu = G.mu(:, idx);
    
    Gs.tag = G.tag;
    Gs.ty = 'm';
    Gs.n = size(mu, 2);
    Gs.d = G.d;
    Gs.mu = mu;
    
    if G.C.n == 1
        Gs.C = G.C;
    else
        Gs.C = pdmat_pick(G.C, idx);
    end

elseif G.ty == 'c'
    
    h = G.h(:, idx);
    
    Gs.tag = G.tag;
    Gs.ty = 'c';
    Gs.n = size(h, 2);
    Gs.d = G.d;
    Gs.h = h;
    
    if G.J.n == 1
        Gs.J = G.J;
    else
        Gs.J = pdmat_pick(G.J, idx);
    end
    
end

