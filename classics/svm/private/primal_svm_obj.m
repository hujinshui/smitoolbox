function [v, g, H] = primal_svm_obj(sol, Q, F, y, c, lof, rb, op)
% The uniform objective function of primal linear/kernel svm training
%
%   The objective is given by  
%
%   minimize    (1/2) w' * Q * w + (rb/2) * b^2 
%             + sum_j c_j L(y_j * (w' * f_j + b)) 
%           
%       w.r.t. (w, b)
%
%   For linear SVM:  Q = I (input 1), and F = X
%   For kernel SVM:  Q = K, and F = K
%
%   [v, g, H] = primal_svm(sol, Q, F, y, c, lof, rb);
%   [v, dx] = primal_svm(sol, Q, F, y, c, lof, rb, 'd'); 
%
%   Inputs:
%       sol:    the solution at which the objective is to be evaluated
%               sol(1:d) - w
%               sol(d+1) - b
%
%       Q:      The quadratic matrix on w [scalar or d x d]
%       F:      The feature matrix [d x n], f_j is given by F(:,j).
%               when F is Q, it is input as empty to note this fact.
%       y:      The label vector [1 x n]
%       c:      The loss coefficients [scalar or 1 x n]
%       lof:    The loss function handle
%       rb:     The regularization coefficient of b [1 x n]
%       dx:     The Newton direction vector: -(H \ g),
%               which is directly computed in an efficient way
%
%       Note for kernel matrix, d = n here.
%
%   Outputs:
%       v:  the objective value [scalar]
%       g:  the gradient vector [(d+1) x 1]
%       H:  the Hessian matrix [(d+1) x (d+1)]
%

% Created by Dahua Lin, on Apr 21, 2011
%

%% parse input

if nargin >= 8 && op(1) == 'd'
    odir = true;
else
    odir = false;
end
    

%% main

% extract solution

if isempty(F)
    d = numel(y);
else
    d = size(F,1);
end
w = sol(1:d);
b = sol(d+1);


% prediction and derivative evaluation

if isempty(F)
    z = y .* (w' * Q + b);
else
    z = y .* (w' * F + b);
end

if nargout <= 1
    Lv = lof(z);
elseif nargout == 2
    if odir
        [Lv, Lg, Lh] = lof(z);
    else
        [Lv, Lg] = lof(z);
    end
else 
    [Lv, Lg, Lh] = lof(z);
end


% evaluate objective

v = eval_objv(Q, rb, w, b, c, Lv);


if odir
    
    % evaluate Newton direction
    if nargout >= 2        
        tg = (c .* Lg) .* y;
        th = c .* Lh;            
        g = eval_dir(Q, F, rb, w, b, tg, th);    
    end
    
else

    % evaluate gradient
    if nargout >= 2
        tg = (c .* Lg) .* y;
        g = eval_grad(Q, F, rb, w, b, tg);
    end
    
    % evaluate Hessian
    
    if nargout >= 3
        th = c .* Lh;
        H = eval_Hess(Q, F, rb, th);
    end
end


%% Sub functions


function v = eval_objv(Q, rb, w, b, c, Lv)

if isscalar(Q)
    vw = (w' * w) * (Q / 2);
else
    vw = (w' * Q * w) / 2;
end
vb = b^2 * (rb / 2);


if isscalar(c)
    v = vw + vb + c * sum(Lv);
else
    v = vw + vb + c * Lv';
end



function g = eval_grad(Q, F, rb, w, b, tg)

gb = rb * b + sum(tg);

if isempty(F)
    g = [Q * (w + tg'); gb];
else
    if isequal(Q, 1)
        Qw = w;
    else
        Qw = Q * w;
    end
    g = [Qw + F * tg'; gb];
end


function H = eval_Hess(Q, F, rb, th)

if isempty(F)
    F = Q;
end
d = size(F, 1);

Hw = F * bsxfun(@times, F, th)';
if isscalar(Q)
    dind = 1 + (0:d-1) * (d+1);
    Hw(dind) = Hw(dind) + Q;
else
    Hw = Hw + Q;
end
Hw = (Hw + Hw') / 2;

Hwb = F * th';

H = zeros(d+1, d+1);

H(1:d, 1:d) = Hw;
H(1:d, d+1) = Hwb;
H(d+1, 1:d) = Hwb';
H(d+1, d+1) = sum(th) + rb;


function dx = eval_dir(Q, F, rb, w, b, tg, th)


if isempty(F) % F == Q
    
    n = size(Q, 1);
    
    ii = [w + tg', th'];
    gw_and_f = Q * ii;  
    gw = gw_and_f(:,1);
    f = gw_and_f(:,2);
        
    do_fast = nnz(th) < 0.7 * n;
    
    if do_fast
        
        % compute M0
        hnz = find(th ~= 0); 
        hz = find(th == 0);
        dhnz = th(hnz); 
        dhnz = dhnz(:);
        M0 = bsxfun(@times, dhnz, Q(hnz, hnz));
        m = size(M0, 1);
        dind = 1 + (0:m-1) * (m+1);
        M0(dind) = M0(dind) + 1;
        
        % compute IA_gw and u
        pre_rr = ii(hnz, :);
        pre_rr(:, 1) = pre_rr(:, 1) - bsxfun(@times, dhnz, Q(hnz, hz) * ii(hz, 1));        
        rr = M0 \ pre_rr;
        IA_gw = ii(:, 1);
        IA_gw(hnz) = rr(:,1);
        u = ii(:, 2);
        u(hnz) = rr(:,2);
        
    else
        
        % compute M
        M = bsxfun(@times, th(:), Q);
        dind = 1 + (0:n-1) * (n+1);
        M(dind) = M(dind) + 1;
        
        % compute IA_gw and u        
        rr = M \ ii;
        IA_gw = rr(:, 1);
        u = rr(:, 2);        
    end

else  % F != Q
    
    d = size(F, 1);
    A = F * bsxfun(@times, F, th)';
    if isscalar(Q)
        dind = 1 + (0:d-1) * (d+1);
        A(dind) = A(dind) + Q;
    else
        A = A + Q;
    end
    A = (A + A') * 0.5;
    
    f = F * th';
    u = A \ f;
    
    gw = Q * w + F * tg';
    IA_gw = A \ gw;
end

gb = rb * b + sum(tg);
hb = rb + sum(th);
s = 1 ./ (hb - f' * u);

db = u' * gw - gb;
dx = [- (IA_gw + (s * db) * u); s * db];

