function R = test_galbp()
% A script to test the correctness of implementation of galbp
%

% Created by Dahua Lin, on Nov 3, 2010
%

%% configure graph

n = 8;
edges = [1 2; 1 8; 2 3; 2 7; 3 4; 3 5; 3 6];
m = size(edges, 1);
w = rand(m, 1);
h = randn(n, 1);

W = zeros(n, n);
W(sub2ind(size(W), edges(:,1), edges(:,2))) = w;
W(sub2ind(size(W), edges(:,2), edges(:,1))) = w;

J = laplacemat(W, rand(n,1));

%% prepare ground-truth

C = inv(J);
v0 = diag(C);
u0 = J \ h;

%% do LBP

bp = galbp(J);

bp.infer_Js('tol', 1e-12);
r_h = bp.infer_hs(h, 'tol', 1e-12);

v_bp = 1 ./ bp.r_Js;
u_bp = r_h ./ bp.r_Js;

%% output

R.J = J;
R.v0 = v0.';
R.u0 = u0.';
R.bp = bp;
R.v_bp = v_bp.';
R.u_bp = u_bp.';
R.v_diff = R.v0 - R.v_bp;
R.u_diff = R.u0 - R.u_bp;

