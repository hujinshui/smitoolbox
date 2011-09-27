function s = dpmm_merge_samples(amdl, ss)
%DPMM_MERGE_SAMPLES Merges multiple samples from the same chain into one
%
%   s = dpmm_merge_samples(amdl, ss);
%   s = dpmm_merge_samples(amdl, ss, rthres);
%
%       merges multiple DPMM samples drawn from a MCMC chain (those
%       produced by the method make_output of dpmm class) into a single
%       sample by voting.
%
%       Input arguments:
%       - amdl:     the underlying atom model
%
%       - ss:       the sequence of samples
%
%       - rthres:   the threshold of ratio. If the total counts of an atom
%                   takes a ratio less than rthres, the atom will be
%                   discarded, and the corresponding observations will be
%                   assigned the labels of the most fitted atom in the
%                   remaining set.
%
%

% Created by Dahua Lin, on Sep 26, 2011.
%

%% verify input arguments

if ~isa(amdl, 'atom_model_base')
    error('dpmm_merge_samples:invalidarg', ...
        'amdl should be an instance of a class derived from atom_model_base.');
end

if ~isstruct(ss) && isvector(ss)
    error('dpmm_merge_samples:invalidarg', 'ss should be a struct vector.');
end



%% main

max_aid = max([ss.max_atom_id]);




