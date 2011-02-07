function [T, sigma0] = test_gatrwa()
% A script to test the correctness of gatrwa
%

% Created by Dahua Lin, on Feb 6, 2011
%

% construct chain graph

n = 300;
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
ca_solve(T, 50000, 1e-9);
fprintf('CA solve: rho diff = %g  sigma diff = %g\n', ...
    norm(rho0 - T.rho), norm(sigma0 - T.sigma));

% F = T.objfunc();
% T.initialize();
% x0 = [T.sigma; T.rho];
% 
% options = optimset('GradObj', 'on', 'LargeScale', 'on',...
%     'MaxIter', 200, 'TolX', 1e-8, 'TolFun', 1e-10, 'Display', 'iter');
% x = fminunc(F, x0, options);
% 
% sigma = x(1:T.nv);
% rho = x(T.nv+1 : end);
% fprintf('Fmin solve: rho diff = %g  sigma diff = %g\n', ...
%     norm(rho0 - rho), norm(sigma0 - sigma));




function ca_solve(T, maxiter, tol)

i = 0;
objv = T.eval_objv();
converged = false;

while ~converged && i < maxiter
            
    i = i + 1;
    
    prev = objv;
     
    T.update_sigma();
    objv = T.eval_objv();
    ch1 = objv - prev;      
    
    prev = objv;
    T.update_rho();
    objv = T.eval_objv();
    ch2 = objv - prev;
    
    ch = ch1 + ch2;        
    fprintf('Iter %d: objv = %f:  ch = %g\n', i, objv, ch);
    
    if abs(ch) < tol
        converged = true;
    end
end

disp(' ');



