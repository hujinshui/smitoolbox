function tf = is_gaussd(G)
% Tests whether the input argument is a gaussd struct.
%
%   tf = is_gaussd(G);
%

% Created by Dahua Lin, on Dec 5, 2011
%

tf = isstruct(G) && isscalar(G) && isfield(G, 'tag') && ...
    strcmp(G.tag, 'gaussd');
