function v = topiclda_ppx(Beta, C, Gamma)
%TOPICLDA_PPX Evaluation of Perplexity for Topic LDA
%
%   v = topiclda_ppx(Beta, C, Gamma);
%
%       Evaluates the perplexity over a corpus for a Topic LDA model.
%
%       Input arguments:
%       - Beta:         The word-distributions of topics
%       - C:            The word count table of the corpus
%       - Gamma:        The per-document posterior Dirichlet param
%                       (which can be solved using topiclda_varinfer)
%
%       Output arguments:
%       - v:            The perplexity over the given document.
%

%% main

Q = exp(psi(Gamma));
Q = bsxfun(@times, Q, 1 ./ sum(Q, 1));
v = topic_ppx(Beta, C, Q);

