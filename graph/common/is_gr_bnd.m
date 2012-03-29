function tf = is_gr_bnd(g)
% Tests whether G is a degree-bounded graph struct
%
%   tf = is_gr_bnd(g);
%       returns whether g is a valid degree-bounded graph struct.
%

% Created by Dahua Lin, on Nov 30, 2011
%

%% main

tf = isstruct(g) && isscalar(g) && isfield(g, 'tag') && strcmp(g.tag, 'gr-bnd');