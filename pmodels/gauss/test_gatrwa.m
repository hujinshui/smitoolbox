function [T, sigma0] = test_gatrwa()
% A script to test the correctness of gatrwa
%

% Created by Dahua Lin, on Feb 6, 2011
%

% construct chain graph

n = 50;
a = zeros(n, 1);
a(1) = 5;
a(n) = 5;
nw = 1.6;

[s,t,w] = gridgraph1d(n, nw);
g = gr_adjlist.from_edges('u', n, s, t, w);
J = laplacemat(g, a);
m = g.ne;

% prepare ground truth

C0 = inv(full(J));
sigma0 = sqrt(diag(C0));

cc = bsxfun(@times, bsxfun(@times, C0, 1 ./ sigma0), 1 ./ sigma0');

rho0 = cc(sub2ind([n n], g.es+1, g.et+1));
rho0 = rho0(1:m);

% initialize model and solve

T = gatrwa(J);
T.set_eprob(1);
T.initialize();

my_solve(T, 10000, 1e-10);

% rho_compare = [rho0, T.rho];
% display(rho_compare);
display(norm(rho0 - T.rho));

% sig_compare = [sigma0, T.sigma];
% display(sig_compare);
display(norm(sigma0 - T.sigma));


function my_solve(T, maxiter, tol)

i = 0;
objv = T.eval_objv();
converged = false;

while ~converged && i < maxiter
    
    fprintf('Iter %d:\n', i);
        
    i = i + 1;
    
    prev = objv;
     
    T.update_sigma();
    objv = T.eval_objv();
    ch1 = objv - prev;
      
    % fprintf('\tupdate sigma: ch = %g\n', ch1);
    
    prev = objv;
    T.update_rho();
    objv = T.eval_objv();
    ch2 = objv - prev;
    % fprintf('\tupdate rho: ch = %g\n', ch2);
    
    ch = ch1 + ch2;        
    fprintf('\tall: objv = %f:  ch = %g\n', objv, ch);
    
    if abs(ch) < tol
        converged = true;
    end
    
    disp(' ');
end

