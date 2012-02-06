function [T, W] = topic_word_count(Beta, C, w, P)
%TOPIC_WORD_COUNT Topic/word counting over corpus
%
%   [T, W] = TOPIC_WORD_COUNT(Beta, C);
%   [T, W] = TOPIC_WORD_COUNT(Beta, C, w);
%   [T, W] = TOPIC_WORD_COUNT(Beta, C, w, P);
%
%       This function measures per-document topic weights and per-topic
%       word weights over a given corpus.
%
%       Specifically, given a corpus with n documents.
%
%       (1) The weight of the k-th topic for the i-th document is defined 
%           to be
%       
%               T(k, i) := sum_v c_i(v) Q_i(k|v).
%
%           Here, c_i(v) is the number of times the word v appears in the 
%           i-th document, and Q_i(k|v), defined below, is the posterior
%           probability that this word belongs to topic k.
%
%               Q_i(k|v) is proportional to p_i(k) * Beta(v, k).
%       
%           such that sum_v Q_i(k|v) = 1 for each v.
%
%
%       (2) The weight of the word v for the k-th topic is defined to be
%
%               W(v, k) := sum_i w_i c_i(v) Q_i(k, v).
%
%       
%       Suppose there are K topics, V distinct words in the vocabulary,
%       and n documents.
%
%       Input arguments:
%       - C:        The word-count matrix of size V x n, where C(v, i)
%                   is the count of word v in the i-th document.
%
%       - w:        The document-weights: a vector of length n.
%                   It can also be empty, indicating that weights are 1. 
%
%       - Beta:     The word-distributions of topics. 
%                   Beta(v, k) is the probability that the word v is
%                   emitted from the k-th topic.
%
%       - P:        The document-specific prior weights of topics.
%                   P(k, i) = p_i(k) is the prior weight of the k-th topic
%                   for the i-th document.
%
%                   P can be a K x n matrix, a K x 1 vector (the same prior
%                   across all documents), or empty (uniform prior).
%
%       Output arguments:
%       - T:        The per-document topic weights, a K x n matrix.
%
%       - W:        The per-topic word weights, a V x K matrix.
%

% Created by Dahua Lin, on Feb 4, 2012
%

%% verify input arguments

if ~(isfloat(Beta) && isreal(Beta) && ndims(Beta) == 2)
    error('topic_word_count:invalidarg', 'Beta should be a real matrix.');
end
[V, K] = size(Beta);

if ~(isfloat(C) && isreal(C) && ndims(C) == 2 && size(C, 1) == V)
    error('topic_word_count:invalidarg', ...
        'C should be a real matrix with V rows.');
end
n = size(C, 2);

if nargin < 3 || isempty(w)
    w = ones(n, 1);
else
    if (isfloat(w) && isreal(w) && isvector(w) && numel(w) == n)
        error('topic_word_count:invalidarg', ...
            'w should be a real vector of length n.');
    end
    if size(w, 2) > 1; w = w.'; end
end

if nargin < 4 || isempty(P)
    P = [];
else
    if ~(isfloat(P) && isreal(P) && ndims(P) == 2 && size(P,1) == K && ...
            (size(P,2) == 1 || size(P,2) == n) )
        error('topic_word_count:invalidarg', ...
            'P should be a matrix of size K x 1 or K x n.');
    end
end

%% main

if isempty(P) || size(P, 2) == 1
    
    Q = Beta.';
    if ~isempty(P)
        Q = bsxfun(@times, Q, P);
    end        
    Q = bsxfun(@times, Q, 1./sum(Q,1));  % Q: K x V
    
    T = Q * C;  % T: K x n
    h = C * w;  % h: V x 1
    W = bsxfun(@times, Q, h.').';
        
else
    
    T = zeros(K, n);
    W = zeros(V, K);
    B = Beta.';     % B: K x V
    
    if ~issparse(C)
        for i = 1 : n            
            c = C(:,i);
            Qi = bsxfun(@times, B, P(:,i));
            Qi = bsxfun(@times, Qi, 1./sum(Qi,1));
            
            T(:,i) = Qi * c;
            W = W + w(i) * bsxfun(@times, Qi, c.').';
        end
    else
        for i = 1 : n
            c = C(:,i);
            [v, ~, c] = find(c);
            Qi = full(bsxfun(@times, B(:,v), P(:,i)));
            Qi = bsxfun(@times, Qi, 1./sum(Qi, 1));
            T(:,i) = Qi * c;
            W(v, :) = W(v, :) + w(i) * bsxfun(@times, Qi, c.').';
        end        
    end                
end
    

