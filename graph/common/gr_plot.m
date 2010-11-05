function [h, coords] = gr_plot(g, coords, varargin)
% Plot a graph
%
%   h = gr_plot(g, 2, ...);
%   h = gr_plot(g, 3, ...);
%   h = gr_plot(g, coords, ...);
%       plots the graph g. 
%
%       In input, g is a struct in form of gr_edgelist (including
%       gr_adjlist), or an affinity matrix. 
%
%       coords is an n x 2 or n x 3 array that gives the coordinates of 
%       the graph (in 2D or 3D space). If coord is 2 or 3, then the
%       function computes a set of 2D or 3D coordinates by spectral
%       embedding. 
%
%       One can additionally specify plotting options.
%
%       This function outputs the handle to the plotted lines.
%
%   [h, coords] = gr_plot( ... );
%       additionally returns the coordinates of the vertices.
%


% Created by Dahua Lin, on Nov 5, 2012
%


%% verify input arguments

g = gr_edgelist(g);

if isscalar(coords)
    d = coords;
    if ~(isequal(d, 2) || isequal(d, 3))
        error('gr_plot:invalidarg', 'The 2nd argument is invalid.');
    end
    coords = [];
else
    if ~(isnumeric(coords) && ndims(coords) == 2 && (size(coords,2) == 2 || size(coords,2) == 3))
        error('gr_plot:invalidarg', ...
            'coords should be an n x 2 or n x 3 numeric matrix.');
    end
end

if nargin < 3
    pargs = {};
else
    pargs = varargin;
end

%% main

% generate coordinates

if isempty(coords)
    coords = gsembed(g, d, 1);
else
    d = size(coords, 2);
end

% generate plot

s = g.s + 1;
t = g.t + 1;

if d == 2
    x = gen_coord_seq(coords(:, 1), s, t);
    y = gen_coord_seq(coords(:, 2), s, t);
    
    h = plot(x, y, pargs{:});
else
    x = gen_coord_seq(coords(:, 1), s, t);
    y = gen_coord_seq(coords(:, 2), s, t);
    z = gen_coord_seq(coords(:, 3), s, t);
    
    h = plot3(x, y, z, pargs{:});
end


%% sub functions        

function xs = gen_coord_seq(x, s, t)

n = length(s);
xs = zeros(n*3, 1);
xs(1:3:end) = x(s);
xs(2:3:end) = x(t);
xs(3:3:end) = nan;


