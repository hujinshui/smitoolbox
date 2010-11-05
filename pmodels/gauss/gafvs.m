function [fvs, scores] = gafvs(cg, nmax)
% Select an (approximate) set of feedback vertices for Gaussian inference
%
%   fvs = gafvs(cg, nmax);
%       selects a set of feedback vertices for Gaussian model inference.
%
%       In input, cg is the partial correlation graph, which can be 
%       constructed by the function parcorrgraph.
%       nmax is a the maximum number of feedback vertices to select.
%
%       The feedback vertices are returned in order of being selected.
%
%   [fvs, scores] = gafvs(cg, nmax);
%       additionally returns the scores when each node is selected.
%

% Created by Dahua Lin, on Nov 5, 2010
%

%% verify input arguments

if ~(isstruct(cg) && isfield(cg, 'tag') && strcmp(cg.tag, 'gr_adjlist') ...
        && strcmp(cg.dty, 'u') && ~isempty(cg.w))
    error('gafvs:invalidarg', 'cg should be an undirected gr_adjlist struct.');
end

if ~(isnumeric(nmax) && isscalar(nmax) && nmax >= 1)
    error('gafvs:invalidarg', 'nmax should be a positive integer scalar.');
end
nmax = double(nmax);

%% main
    
if nargin < 2
    fvs = gafvs_cimp(cg, nmax);
else
    [fvs, scores] = gafvs_cimp(cg, nmax);
end

