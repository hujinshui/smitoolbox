function test_agmrf()
% A simple script to test the implementation of gmrf and agmrf
%

% Created by Dahua Lin, on Feb 12, 2011
%

%% configure model

n = 15;
ker = [1.2, 0.8, 0.4];

[s, t, w] = gridgraph1d(n, ker);
g = gr_adjlist.from_edges('u', n, s, t, w);

a = rand(n, 1);

agm = agmrf(g, a);
gm = to_gmrf(agm);

x = randn(n, 1);
y = randn(n, 1);


%% compute information matrix

disp('Information matrix:');

J0 = full(laplacemat(g, a));
Ja = information_matrix(agm);
Jg = information_matrix(gm);

fprintf('isequal(J0, Ja, Jg) = %d\n', isequal(J0, Ja, Jg));

%% compute energy

disp('Energy evaluation:');

int_e0 = (x' * laplacemat(g) * x) / 2;
int_ea = agm.int_energy(x);
fprintf('|int_ea - int_e0| = %g\n', abs(int_ea - int_e0));

ext_e0 = a' * (x - y).^2 / 2;
ext_ea = agm.ext_energy(x, y);
fprintf('|ext_ea - ext_e0| = %g\n', abs(ext_ea - ext_e0));

tot_e0 = int_e0 + ext_e0;
tot_ea = agm.energy(x, y);
fprintf('|tot_ea - tot_e0| = %g\n', abs(tot_ea - tot_e0));


%% gmrf expected energy

disp('gmrf energy evaluation:');

C = randn(n); C = C * C';

vs = diag(C);
es = g.es(1:g.ne) + 1;
et = g.et(1:g.ne) + 1;

cs = C(sub2ind([n n], es, et));

sigma = sqrt(vs);
rho = cs ./ (sigma(es) .* sigma(et));

ge0 = trace(J0 * C) / 2;
ge = gm.energy(vs, cs);
ge_sr = gm.energy_sr(sigma, rho);

fprintf('|ge - ge0| = %g\n', abs(ge - ge0));
fprintf('|ge_sr - ge0| = %g\n', abs(ge_sr - ge0));


%% expected energy evaluation

disp('Expected energy evaluation:');

h = a .* y;
e0 = ge0 + (x' * J0 * x) / 2 + (h' * y) / 2 - h' * x;
eg = ge + gm.qterm(x) + (h' * y) / 2 - h' * x;
ea = agm.energy_ep(x, vs, cs, y);

fprintf('|eg - e0| = %g\n', abs(eg - e0));
fprintf('|ea - e0| = %g\n', abs(ea - e0));


