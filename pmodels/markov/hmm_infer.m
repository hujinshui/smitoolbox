function [Q, H] = hmm_infer(Pi, T, LL)
%HMM_INFER HMM E-step inference
%
%   Q = HMM_INFER(Pi, T, LL);
%
%       Infers the (posterior) marginal distributions of all states in 
%       an Hidden Markov Model.
%
%       Suppose there are K distinct state values, and n steps along 
%       the chain.
%
%       Input arguments:
%       - Pi:       The initial distribution of states [K x 1]
%       - T:        The transition probability matrix [K x K]
%       - LL:       The matrix of log-likelihood values [K x n]
%
%       In the output, Q is a matrix of K x n, where Q(k, i) is the
%       posterior probability of the state at i-th step being k,  
%       given all observations.
%
%   [Q, H] = HMM_INFER(Pi, T, LL);       
%
%       Also computes H, the accumulated transition counts.
%
%       Specifically, H is an K x K matrix, with
%
%       H(u, v) = sum_{i=2}^n pr(x_{n-1}=u, x_n=v | all observations)
%
%   Remarks
%   -------
%       - This function actually performs the E-step in E-M estimation 
%         of an Hidden markov model.
%

% Created by Dahua Lin, in Feb 1, 2012
%

%% main

% run forward-backward recursion

[A, Lc] = hmm_forward(Pi, T, LL);   % hmm_forward will verify all inputs
B = hmm_backward(T, LL, Lc);

% extract outputs

Q = A .* B;

% calculate H

if nargout < 2; return; end

E = exp(bsxfun(@minus, LL(:,2:end), Lc(2:end)));
H = (A(:,1:end-1) * (E .* B(:,2:end))') .* T;

