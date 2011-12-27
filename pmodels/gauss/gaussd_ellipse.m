function gaussd_ellipse(G, r, n, varargin)
% Draws an Ellipse to represent a Gaussian model
%
%   gaussd_ellipse(G, r, n, ...);
%
%   Input arguments:
%       - G:        A gaussd struct to represent the Gaussian model(s)
%       - r:        the relative radius
%       - n:        the number of points on each ellipse
%
%   Output:
%       - h:        the handle to the drawn lines.
%
%   One can also specify additional options to customize the plot like
%   as using the plot function.
%

% Created by Dahua Lin, on Dec 27, 2011
%

%% verify input

if ~(is_gaussd(G) && G.d == 2)
    error('gaussd_ellipse:invalidarg', ...
        'G should be a gaussd struct with G.d == 2.');
end

if ~(isfloat(r) && isreal(r) && isscalar(r) && r > 0)
    error('gaussd_ellipse:invalidarg', ...
        'r should be a positive real value.');
end

if ~(isnumeric(n) && isscalar(n) && n == fix(n) && n >= 8)
    error('gaussd_ellipse:invalidarg', ...
        'n should be a positive integer with n >= 8.');
end

%% main

t = linspace(0, 2*pi, n);
Z = [cos(t); sin(t)];
if r ~= 1
    Z = Z * r;
end


if G.n == 1
    X = gaussd_sample(G, [], Z);
    plot(X(1,:), X(2,:), varargin{:});
else
    for i = 1 : G.n
        g = gaussd_sub(G, i);
        X = gaussd_sample(g, [], Z);
        plot(X(1,:), X(2,:), varargin{:});
        hold on;
    end
    hold off;
end

