function [T, sigma0] = test_gatrwa(siz)
% A script to test the correctness of gatrwa
%

% Created by Dahua Lin, on Feb 6, 2011
%


if numel(siz) == 1
    [T, sigma0] = test1d(siz);
elseif numel(siz) == 2
    [T, sigma0] = test2d(siz);
end



function [T, sigma0] = test1d(n)

% construct chain graph

a = zeros(n, 1);
a(1) = 5;
a(n) = 5;
a = a + 0.01 * rand(size(a));
nw = 1.6;

[s,t,w] = gridgraph1d(n, nw);
g = gr_adjlist.from_edges('u', n, s, t, w);
J = laplacemat(g, a);

% solve

[T, sigma0] = do_task(J, g, 1);

% visualize

plot(1:T.nv, sigma0, 1:T.nv, T.sigma);


function [T, sigma0] = test2d(siz)

% construct mrf grid

n1 = siz(1);
n2 = siz(2);

a = zeros(n1, n2);
ea = 5;
a(1, :) = ea;
a(n1, :) = ea;
a(:, 1) = ea;
a(:, n2) = ea;
nw = 1.6;

[s,t,w] = gridgraph2d([n1 n2], [0, nw; nw, 0]);
g = gr_adjlist.from_edges('u', n1 * n2, s, t, w);
J = laplacemat(g, a(:));

% solve

[T, sigma0] = do_task(J, g, 0.75);

% visualize

plot(1:T.nv, sigma0, 1:T.nv, T.sigma);



function [T, sigma0, rho0] = do_task(J, g, ep)

% prepare ground truth

n = size(J, 1);
m = g.ne;

C0 = inv(full(J));
sigma0 = sqrt(diag(C0));

cc = bsxfun(@times, bsxfun(@times, C0, 1 ./ sigma0), 1 ./ sigma0');

rho0 = cc(sub2ind([n n], g.es+1, g.et+1));
rho0 = rho0(1:m);

% initialize model and solve

T = gatrwa(J);
T.set_eprob(ep);

T.initialize();
T.solve('maxiter', 50000, 'display', true);
fprintf('Compare with groud-truth: rho diff = %g  sigma diff = %g\n', ...
    norm(rho0 - T.rho) / sqrt(numel(rho0)),  ...
    norm(sigma0 - T.sigma) / sqrt(numel(sigma0)));
