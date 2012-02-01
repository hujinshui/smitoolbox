function T = markov_tpm(K, seqs, w, pri_c)
%MARKOV_TPM Estimation of Markov Transition Probability Matrix
%
%   T = MARKOV_TPM(K, seqs);
%   T = MARKOV_TPM(K, seqs, w);
%   T = MARKOV_TPM(K, seqs, [], pri_c);
%   T = MARKOV_TPM(K, seqs, w, pri_c);
%
%       Estimates the transition probability from observed sequences
%
%       Input arguments:
%       - K:        The number of different states
%
%       - seqs:     The sequences, which can be in either of the following
%                   form. 
%
%                   - a row vector of states
%
%                   - a K x n matrix, seqs(k, i) gives the probability that
%                     the i-th step is in the k-th state.
%                     Each column of this matrix needs to sum to 1.
%
%                   - a cell array. The element in each cell is in either
%                     of the above forms.
%
%       - w:        The weights of the sequences. Suppose there are
%                   m sequences, then w should be a vector of length m.
%                   
%       - pri_c:    The prior count. (If pri_c is given, it assumes a 
%                   Dirichlet prior with alpha = pri_c + 1, and MAP
%                   estimation is to be performed.)
%
%                   pri_c can be either a scalar or an K x K matrix.
%

% Created by Dahua Lin, on Jan 31, 2012
%

%% verify input arguments

if ~(isnumeric(K) && isscalar(K) && K >= 2 && K == fix(K))
    error('markov_tpm:invalidarg', 'K must be an integer with K >= 2.');
end

if isnumeric(seqs)
    m = 1;
else
    m = numel(seqs);
    if m == 1
        seqs = seqs{1};
    end
end

if nargin >= 3 && ~isempty(w)
    if ~(isfloat(w) && isreal(w) && numel(w) == m)
        error('markov_tpm:invalidarg', ...
            'w should be a reall array with m elements.');
    end
else
    w = [];
end

if nargin >= 4
    if ~(isfloat(pri_c) && ...
            ((isscalar(pri_c) && pri_c >= 0) || isequal(size(pri_c), [K K]))  )
        error('markov_tpm:invalidarg', ...
            'pri_c should be a non-negative real scalar or a K x K real matrix.');
    end
else
    pri_c = 0;
end

%% main

% summarize observations

if m == 1
    sH = count_ts(K, seqs);
    if ~isempty(w) && w ~= 1
        sH = sH * w;
    end
else
    sH = zeros(K, K);    
    if isempty(w)                
        for i = 1 : m
            H = count_ts(K, seqs{i});
            sH = sH + H;
        end
    else        
        for i = 1 : m
            H = count_ts(K, seqs{i});
            sH = sH + H * w(i);
        end
    end
end

% solve estimation

if ~isequal(pri_c, 0)
    sH = sH + pri_c;
end

T = bsxfun(@times, sH, 1 ./ sum(sH, 2));


%% counting functions

function H = count_ts(K, s)

if ~(isnumeric(s) && isreal(s) && ndims(s) == 2)
    error('markov_tpm:invalidarg', ...
        'Each sequence must be a real vector/matrix.');
end

if size(s, 1) == 1
    I = s(1:end-1);
    J = s(2:end);
    L = I + (J - 1) * K;
        
    H = intcount(K*K, L);
    H = reshape(H, K, K);
    
elseif size(s, 1) == K
    A = s(:, 1:end-1);
    B = s(:, 2:end);
    H = A * B';
    
else
    error('markov_tpm:invalidarg', ...
        'The size of some sequence is invalid.');
end

