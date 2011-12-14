function R = gaussd_mapest(G, dh, dJ)
% Performs MAP estimation w.r.t. Gaussian prior
%
%   R = gaussd_mapest(G, dh, dJ);
%
%       Performs MAP estimation with the Gaussian prior G, and the 
%       updates to the canonical param of the prior.
%       
%       Input arguments:
%       - G:    the Gaussian prior (G.ty == 'c' && G.n == 1)
%       - dh:   the update to the potential vector
%       - dJ:   the update to the precision matrix.
%
%       Output arguments:
%       - R:    the MAP estimation.
%
%       Generally, R is given by (J_pri + dJ) \ (h_pri + dh).
%       

% Created by Dahua Lin, on Dec 14, 2011
%

%% main

[h, J] = gaussd_conjupdate(G, dh, dJ);
R = pdmat_lsolve(J, h);      

