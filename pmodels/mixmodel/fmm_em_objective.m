function v = fmm_em_objective(pri, LL, w, c0, s)
%FMM_EM_OBJECTIVE Evaluate the objective of an FMM-EM solution
%
%   v = FMM_EM_OBJECTIVE(pri, LL, w, c0, s);
%
%       Evaluates the objective function value for EM-estimation of
%       a finite mixture model.
%

% Created by Dahua Lin, on Feb 3, 2012
%

%% main

if ~isempty(w)
    w = w(:);
end

% log-pri: Pi

log_pi = log(s.Pi);

if isequal(c0, 0)
    lpri_pi = 0;
else
    lpri_pi = c0 * sum(log_pi);
end

% log-pri: params

if isempty(pri)
    lpri_t = 0;
else
    lpri_t = pri.logpdf(s.params);
    lpri_t = sum(lpri_t);
end

% log-lik: labels (Z)

Q = s.Q;
llik_z = sum_w(log_pi' * Q, w);

% log-lik: observations

llik_x = sum_w(sum(Q .* LL, 1), w);

% entropy

ent = sum_w(ddentropy(Q), w);

% overall sum

v = lpri_pi + lpri_t + llik_z + llik_x + ent;


%% auxiliary function

function v = sum_w(x, w)

if isempty(w)
    v = sum(x);
else
    v = x * w;
end





