function tf = is_gr(g)
% Tests whether the input argument is a valid graph struct
%
%   tf = is_gr(g);
%
%       Tests whether the input argument g is a valid graph struct 
%       produced by make_gr.
%

% Created by Dahua Lin, on Oct 28, 2011
%

%% main

tf = isstruct(g) && isscalar(g) && isfield(g, 'tag') && strcmp(g.tag, 'gr');
