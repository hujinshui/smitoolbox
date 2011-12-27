function R = gaussd_mapest(G, S)
% Performs MAP estimation w.r.t. Gaussian prior
%
%   R = gaussd_mapest(G, S);
%
%       Performs MAP estimation with the Gaussian prior G, and the 
%       updates to the canonical param of the prior.
%       
%       Input arguments:
%       - G:    the Gaussian prior (G.ty == 'c' && G.n == 1)
%       - S:    the gaussd struct that captures the updates derived
%               from the observations.
%
%       Output arguments:
%       - R:    the MAP estimation.
%
%       Generally, R is given by (J_pri + dJ) \ (h_pri + dh).
%       

% Created by Dahua Lin, on Dec 14, 2011
%

%% main

Gp = gaussd_conjupdate(G, S);
R = pdmat_lsolve(Gp.J, Gp.h);      

