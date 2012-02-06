function [objv, s] = topiclda_varinfer_objv(Beta, alpha, Gamma, sumPhi, entPhi)
%TOPICLDA_VARINFER_OBJV Objective function of LDA Variational Inference
%
%   [objv, s] = TOPICLDA_VARINFER_OBJV(Beta, alpha, Gamma, sumPhi, entPhi);
%

%% main

K = size(Beta, 2);

% expected log-lik of theta (w.r.t. gamma)

A = sum(psi(Gamma), 1) - psi(sum(Gamma, 1));  
ell_theta = gammaln(K * alpha) - K * gammaln(alpha) + (alpha - 1) * A;
ell_theta = sum(ell_theta);

% expected log-lik of z (w.r.t gamma and 


