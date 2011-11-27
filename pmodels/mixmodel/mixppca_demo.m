function mixppca_demo(n, K)
% A program to demonstrate the use of mixture of PPCA
%
%   mixppca_demo(n, K);
%
%       Here, n is the number of sample points, K is the number of 
%       mixture components to be used
%

% Created by Dahua Lin, on Nov 6, 2011
%

%% generate data

sigma0 = 0.1;
r0 = 2;
rgn = [-3, 3, -3, 3];

theta = rand(1, n) * (2 * pi);
r = r0 + randn(1, n) * sigma0;

X = [r .* cos(theta); r .* sin(theta)];

%% Do estimation

gm = ppca_gm(2, 1);
prg = fmm_std(gm, [], K);

opts = varinfer_options([], ...
    'maxiters', 500, ...
    'display', 'eval', ...
    'tol', 1e-9);

R = smi_varinfer(prg, X, [], opts);
models = R.sol.params;


%% visualize

% models

figure;
plot(X(1,:), X(2,:), 'b.', 'MarkerSize', 3);

axis(rgn);
axis equal;

for k = 1 : K
    G = models{k}.to_gauss('mp');
    hold on;
    G.plot_ellipse(2, 'r');
end

% likelihood map

tx = linspace(rgn(1), rgn(2), 512);
ty = linspace(rgn(3), rgn(4), 512);
[xx, yy] = meshgrid(tx, ty);

Xg = [xx(:)'; yy(:)'];

P = zeros(K, size(Xg, 2));
for k = 1 : K
    P(k, :) = models{k}.pdf(Xg);
end
P = R.sol.Pi(:)' * P;

P = reshape(P, size(xx));
figure;
imagesc(tx, ty, P);
axis equal;


