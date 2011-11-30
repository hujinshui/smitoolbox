function g = gr_nbs(g)
% Add the neighborhood system to a graph
%
%   g = gr_nbs(g);
%
%       This function creates or updates the neighborhood system for g.
%

% Created by Dahua Lin, on Nov 30, 2011
%

if ~is_gr(g)
    error('gr_nbs:invalidarg', 'g should be a graph struct.');
end

[g.o_nbs, g.o_eds, g.o_degs, g.o_os] = gr_nbs_cimp(g);    
g.has_nbs = true;

