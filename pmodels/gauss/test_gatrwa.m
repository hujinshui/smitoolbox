function [GM, sigma0, sigma] = test_gatrwa(siz)
% A script to test the correctness of gatrwa
%

% Created by Dahua Lin, on Feb 6, 2011
%


if numel(siz) == 1
    [GM, sigma0, sigma] = test1d(siz);
elseif numel(siz) == 2
    [GM, sigma0, sigma] = test2d(siz);
end



function [GM, sigma0, sigma] = test1d(n)

% construct chain graph

a = zeros(n, 1);
a(1) = 5;
a(n) = 5;
a = a + 0.01 * rand(size(a));
nw = 1.6;

[s,t,w] = gridgraph1d(n, nw);
g = gr_adjlist.from_edges('u', n, s, t, w);

% solve

[GM, sigma0, sigma] = do_task(g, a, 1);

% visualize

plot(1:GM.nv, sigma0, 1:GM.nv, sigma);


function [GM, sigma0, sigma] = test2d(siz)

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

% solve

[GM, sigma0, sigma] = do_task(g, a(:), 0.75);

% visualize

plot(1:GM.nv, sigma0, 1:GM.nv, sigma);



function [GM, sigma0, sigma] = do_task(g, a, ep)

% prepare ground truth

n = g.nv;

GM = sgmrf.amodel(g, a);
J = infomat(GM);
C0 = inv(full(J));
vs0 = diag(C0);
cs0 = C0(sub2ind([n n], GM.es, GM.et));

% solve

[vs, cs] = gatrwa(GM, ep, 'display', true);
fprintf('Compare with groud-truth: vs diff = %g, cs diff = %g\n', ...
    mean(abs(vs - vs0)), mean(abs(cs - cs0)) );

% output

sigma0 = sqrt(vs0).';
sigma = sqrt(vs).';



