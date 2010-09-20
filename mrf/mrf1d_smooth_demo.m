function mrf1d_smooth_demo(x0, noise_mag, mker)
% A function to demo 1D curve smoothing based on 1D MRF
%
%   mrf1d_smooth_demo(x0, noise_mag, mker);
%
%   Inputs:
%       x0:         the sequence of values (not corruped by noise)
%       noise_mag:  the magnitude of noise
%       mker:       the mrf kernel
%

% Created by Dahua Lin, on Sep 19, 2010
%

%% verify input

if ~(isfloat(x0) && isvector(x0))
    error('mrf1d_smooth_demo:invalidarg', ...
        'x0 should be a numeric vector.');
end

if ~(isfloat(noise_mag) && isscalar(noise_mag))
    error('mrf1d_smooth_demo:invalidarg', ...
        'noise_mag should be a numeric scalar.');
end

if nargin < 3
    mker = 1;
end

%% main

% prepare data

n = length(x0);
x = x0 + randn(size(x0)) * noise_mag;

% construct MRF

W = mrf1d(n, mker);
M = L2mrf(W);

% do smooth

xs = M.smooth(x, 1);

% visualize

figure;
plot(x, 'b.');
hold on; plot(x0, 'g-');
hold on; plot(xs, 'r-');


