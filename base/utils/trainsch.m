function sch = trainsch(n, n0, gr, rstream)
% Makes a training schedule
%
%   sch = trainsch(n, n0, gr);
%   sch = trainsch(n, n0, gr, rstream);
%
%       In the practice of modeling training, it is sometimes more 
%       efficient to train a model in several stages. Starting from
%       a small random subset to train an initial model, and use it
%       as an initial solution for further optimization with large
%       sample set.
%
%       This function is to produce such a training schedule. 
%       
%       Input arguments:
%       - n:        the total number of available training samples
%       - n0:       the size of initial training set
%       - gr:       the grow ratio. For example, if gr = 2, the training
%                   set doubles at each stage. (gr > 1)
%       - rstream:  the random stream for producing the random numbers.
%                   (this can be omitted).
%       
%       Output argument:
%       - sch:      It is a cell array of size K x 1, where K is the 
%                   number of stages. sch{k} is an array of indices
%                   selected for the k-th stage training. 
%                   sch{K} is always [], which indicates that the entire
%                   training set is used for the final training.
%
%                   For each k > 1, sch{k} is a super-set of sch{k-1}.
%
%       Note if n0 <= n, sch will be {[]}, which indicates there is just
%       one-stage training with all samples.
%

%   History
%   -------
%       - Created by Dahua Lin, on April 24, 2011
%

%% verify input arguments

if ~(isnumeric(n) && isscalar(n) && n > 0 && n == fix(n))
    error('trainsch:invalidarg', 'n should be a positive integer scalar.');
end

if ~(isnumeric(n0) && isscalar(n0) && n0 > 0 && n0 == fix(n0))
    error('trainsch:invalidarg', 'n0 should be a positive integer scalar.');
end

if ~(isnumeric(gr) && isscalar(gr) && isreal(gr) && gr > 1)
    error('trainsch:invalidarg', 'gr should be a real value with gr > 1.');
end

if nargin < 4
    rstream = [];
else
    if ~isa(rstream, 'RandStream')
        error('trainsch:invalidarg', 'rstream should be an object of class RandStream.');
    end
end

%% main

if n > n0
    K = ceil(log(n / n0) / log(gr)) + 1;
    ns = ceil(n0 * (gr.^(0:K-1)));
    ds = [ns(1), diff(ns)];
    
    sch = cell(K, 1);
    
    rs = 1 : n;
    
    for k = 1 : K-1
        a = randpick(numel(rs), ds(k), rstream);
        if k == 1
            s = a;
        else
            s = [s, rs(a)]; %#ok<AGROW>
        end
        sch{k} = s;        
        rs(a) = [];
    end    
    
    sch{K} = [];    
else
    sch = {[]};
end
    




