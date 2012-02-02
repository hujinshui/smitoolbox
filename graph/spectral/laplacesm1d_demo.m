function laplacesm1d_demo(x0, sigma, w)
% A function to demo 1D Laplacian smoothing
%
%   laplacesm1d_demo(x0, noise_mag, mker);
%
%   Inputs:
%       x0:         the sequence of values (not corruped by noise)
%       sigma:      the standard deviation of noise
%       w:          the neighboring link weights (size 1 x r)
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

if ~(isfloat(sigma) && isscalar(sigma))
    error('laplacesm1d_demo:invalidarg', ...
        'sigma should be a numeric scalar.');
end

if nargin < 3
    w = 1;
else
    if ~(isfloat(w) && isreal(w) && ndims(w) == 2 && size(w,1) == 1)
        error('laplacemat1d_demo:invalidarg', ...
            'w should be a real row vector.');
    end
    w = w.';
end
r = numel(w);


%% main

% prepare data

n = length(x0);
x = x0 + randn(size(x0)) * sigma;

% construct MRF

[gb, dx] = gr_local(n, r);  
ws = repmat(w(abs(dx)), 1, n);
[Gs, ws] = gr_sym(gb, ws, 'max');

% do smooth

a = 1;
xs = laplacesm(Gs, ws, a, x);

% visualize

figure;
plot(x, 'b.', 'MarkerSize', 3);
hold on; plot(x0, 'g-');
hold on; plot(xs, 'r-');


