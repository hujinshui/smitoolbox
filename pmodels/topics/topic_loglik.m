function LL = topic_loglik(Beta, C, Q)
%TOPIC_LOGLIK Document log-likelihood for Topic model
%
%   LL = TOPIC_LOGLIK(Beta, C, Q);
%       
%       Evaluates the log-likelihood of each document based on a
%       topic model.
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
%       Output arguments:
%       - LL:       The log-likelihood of all documents [1 x n].
%
%   Remarks
%   -------
%       - This is a generic function that can help to evaluate perplexity
%         for different types of topic models.
%

% Created by Dahua Lin, on Feb 18, 2012.
%

%% verify input arguments

if ~(isfloat(Beta) && isreal(Beta) && ismatrix(Beta) && ~issparse(Beta))
    error('topic_loglik:invalidarg', 'Beta should be a non-sparse real matrix.');
end
[V, K] = size(Beta);

if ~(isfloat(C) && isreal(C) && ismatrix(C) && size(C,1) == V)
    error('topic_loglik:invalidarg', 'C should be a real matrix with V rows.');
end
n = size(C, 2);

if ~(isfloat(Q) && isreal(Q) && isequal(size(Q), [K n]) && ~issparse(Q))
    error('topic_loglik:invalidarg', ...
        'Q should be a non-sparse real matrix of size K x n.');
end

%% main

[I, J, wc] = find(C);
I = int32(I) - 1;
J = int32(J) - 1; 

if ~isa(Beta, 'double'); Beta = double(Beta); end
if ~isa(wc, 'double'); wc = double(wc); end
if ~isa(Q, 'double'); Q = double(Q); end

LL = topic_loglik_cimp(Beta, I, J, wc, Q);