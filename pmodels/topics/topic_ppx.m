function v = topic_ppx(Beta, C, Q)
%TOPIC_PPX Corpus Perplexity for Topic model
%
%   v = TOPIC_PPX(Beta, C, Q);
%       
%       Evaluates the perplexity of a topic model on a corpus.
%
%       Suppose there are K topics, n documents and V words in the
%       vocabulary.
%
%       Input arguments:
%       - Beta:     The word distributions of topics [V x K]
%       - C:        The word count table of the corpus [V x n]
%       - Q:        The per-doc posterior topic distribution [K x n]
%       
%
%       This function outputs v, the perplexity over the corpus.
%
%   Remarks
%   -------
%       - This is a generic function that can help to evaluate perplexity
%         for different types of topic models.
%

% Created by Dahua Lin, on Feb 18, 2012.
%

%% main

LL = topic_loglik(Beta, C, Q);
tc = full(sum(C(:)));
v = exp(sum((-1/tc) * LL));


