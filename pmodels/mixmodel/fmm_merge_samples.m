function s = fmm_merge_samples(pri, gm, c0, X, w, ss)
%FMM_MERGE_SAMPLES Merge multiple FMM solutions into a single one
%
%   s = FMM_MERGE_SAMPLES(pri, gm, X, w, ss);
%
%       Creates a single FMM that best represents all input FMM samples.
%
%       Input arguments:
%       - pri:          the prior object
%       - gm:           the generative model
%       - X:            the observations
%       - w:            the observation weights
%       - ss:           the cell array of FMM samples
%
%       Output arguments:
%       - s:            the merged FMM sample.
%

% Created by Dahua Lin, on Feb 3, 2011
%

%% verify input arguments

if ~iscell(ss)
    error('fmm_merge_samples:invalidarg', ...
        'ss should be a cell array of FMM samples.');
end

%% main

K = ss{1}.K;
n = numel(ss);

% get Z through majority voting

Zm = cell(n, 1);
for i = 1 : n
    cs = ss{i};
    if cs.K ~= K
        error('fmm_merge_samples:invalidarg', ...
            'The K values are inconsistent.');
    end    
    Zm{i} = cs.z;
end
Zm = vertcat(Zm{:});
z = mode(Zm, 1);

% re-estimate parameters

params = fmm_est_params(pri, gm, X, w, {K, z});
Pi = ddestimate({K, z}, w, c0);

% make merged sample

s.K = K;
s.Pi = Pi;
s.params = params;
s.z = z;

            