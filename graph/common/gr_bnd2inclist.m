function G = gr_bnd2inclist(G0)
% Convert a degree-bounded graph struct to a incidence list struct
%
%   G = gr_bnd2inclist(G0);
%
%       converts a degree-bounded graph struct G0 to a directed graph
%       struct with neighborhood system.
%
%       Note that if G0 is associated with a weight matrix W, the
%       corresponding weight vector for G is W(G0.nbs > 0);
%

% Created by Dahua Lin, on Nov 30, 2011
%

%% verify input

if ~is_gr_bnd(G0)
    error('gr_bnd2inclist:invalidarg', ...
        'G0 should be a degree-bounded graph struct.');
end

%% main

n = G0.n;
K = G0.K;

s = repmat(1:n, K, 1);
t = G0.nbs;

s = s(t > 0);
t = t(t > 0);

G = make_gr('d', n, s, t, 'nbs');


