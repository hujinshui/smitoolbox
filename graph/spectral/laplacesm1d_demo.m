function laplacesm1d_demo(x0, noise_mag, gker)
% A function to demo 1D Laplacian smoothing
%
%   laplacesm1d_demo(x0, noise_mag, mker);
%
%   Inputs:
%       x0:         the sequence of values (not corruped by noise)
%       noise_mag:  the magnitude of noise
%       gker:       the grid graph kernel
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 19, 2010
%       - Modified by Dahua Lin, on Nov 2, 2010
%       - Modified by Dahua Lin, on Nov 13, 2010
%

%% verify input

if ~(isfloat(x0) && isvector(x0))
    error('laplacesm1d_demo:invalidarg', ...
        'x0 should be a numeric vector.');
end
if size(x0, 2) > 1
    x0 = x0.';
end

if ~(isfloat(noise_mag) && isscalar(noise_mag))
    error('laplacesm1d_demo:invalidarg', ...
        'noise_mag should be a numeric scalar.');
end

if nargin < 3
    gker = 1;
end

%% main

% prepare data

n = length(x0);
x = x0 + randn(size(x0)) * noise_mag;

% construct MRF

[s,t,w] = gridgraph1d(n, gker);
g = gr_adjlist.from_edges('u', n, s, t, w);
a = 1;

% do smooth

xs = laplacesm(g, a, x);

% visualize

figure;
plot(x, 'b.');
hold on; plot(x0, 'g-');
hold on; plot(xs, 'r-');


