function [vs, ds] = gr_flood(G, s, msk)
%GR_FLOOD Breadth-first flooding
%
%   [vs, ds] = GR_FLOOD(G, s);
%   [vs, ds] = GR_FLOOD(G, s, msk);
%
%       Performs a breadth-first traversal from specified sources (with
%       in allowable region)
%
%       Input arguments:
%       - G:        The graph (with neighbor system)
%
%       - s:        The source of vector of sources
%
%       - msk:      The mask of region of interests, which should be
%                   a logical vector of length G.n.
%
%                   One can set msk(boundary) to false in order to prevent
%                   the traversal beyond some boundary.
%
%       Output arguments:
%       - vs:       The vector of vertices listed in the order of visiting.
%                   (not including the sources)
%
%       - ds:       The vector of corresponding distances to sources.
%                   

% Created by Dahua Lin, on Jan 28, 2012
%

%% verify inputs

if ~(is_gr(G) && G.has_nbs)
    error('gr_flood:invalidarg', ...
        'G should be a graph (with neighbor system).');
end
n = G.n;

if ~(isnumeric(s) && isvector(s) && isreal(s) && ~issparse(s))
    error('gr_flood:invalidarg', 's should be a numeric vector.');
end
s = int32(s);

if ~all(s >= 1 & s<= n)
    error('gr_flood:invalidarg', 'Some source vertices are not valid.');
end

if nargin < 3 || isempty(msk)
    msk = [];
else
    if ~(islogical(msk) && numel(msk) == n)
        error('gr_flood:invalidarg', ...
            'msk should be a logical matrix with n elements.');
    end
end

%% main

if nargout <= 1
    vs = gr_flood_cimp(G, s, msk);
else
    [vs, ds] = gr_flood_cimp(G, s, msk);
end




