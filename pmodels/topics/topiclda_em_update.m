function [sol, objVs] = topiclda_em_update(C, w, eta, sol, eiters, objVs)
%TOPICLDA_EM_UPDATE Performs E-M update of a Topic-LDA solution
%
%   [sol, objVs] = TOPICLDA_EM_UPDATE(C, w, eta, sol, eiters);
%
%       Performs E-M update for Topic-LDA estimation.
%
%       Input arguments:
%       - C:        The word-count matrix [V x n]
%       - w:        The document weights [empty or 1 x n]
%       - eta:      The prior-count for word-distributions
%       - sol:      The topiclda solution to be updated
%       - eiters;   The maximum iterations within each E-step
%
%       Output arguments:
%       - sol:      The updated solution struct
%       - objVs:    The itemized objective values
%
%
%   Remarks
%   -------
%       - This function is to help the implementation of topic-LDA,
%         and is not supposed to be directly called by end users.
%

% Created by Dahua Lin, on Feb 18, 2012
%

%% main

alpha = sol.alpha;
Beta = sol.Beta;
Gamma = sol.Gamma;

% E-step (variational inference)

if isempty(Gamma)
    sC = full(sum(C, 1));
    Gamma = bsxfun(@plus, alpha, sC * (1/sol.K));
end

[Gamma, W, Vs] = topiclda_varinfer(Beta, alpha, C, w, Gamma, ...
    eiters, 1e-12);

if isempty(w)
    Vs = sum(Vs, 2);
else
    Vs = Vs * w';
end

objVs.ell_theta = Vs(1);
objVs.ell_z = Vs(2);
objVs.ell_w = Vs(3);
objVs.ent_theta = Vs(4);
objVs.ent_z = Vs(5);

% M-step (parameter estimation)

% update alpha

sg = sum(Gamma, 1);
dird_sta = bsxfun(@minus, psi(Gamma), psi(sg));
alpha = dird_mle(dird_sta, w, alpha, 'input', 'stat');

et = calc_eloglik_theta(alpha, dird_sta);
objVs.ell_theta = sum(et);

% update Beta

Wt = W.';
if eta > 0
    W1 = Wt + eta;
else
    W1 = Wt';
end

Beta = bsxfun(@times, W1, 1 ./ sum(W1, 1));

logB = log(Beta);
objVs.ell_w = sum(sum(Wt .* logB, 1));

if eta > 0
    objVs.lpri_beta = eta * sum(logB(:));
end


% Store updated results

sol.alpha = alpha;
sol.Beta = Beta;
sol.Gamma = Gamma;
sol.W = W;


%% Aux functions

function v = calc_eloglik_theta(alpha, dird_sta)

v = gammaln(sum(alpha)) - sum(gammaln(alpha));
v = v + (alpha - 1)' * dird_sta;

