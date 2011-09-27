function s = dpmm_merge_samples(amdl, obs, H, ss, rthres)
%DPMM_MERGE_SAMPLES Merges multiple samples from the same chain into one
%
%   s = dpmm_merge_samples(amdl, obs, H, ss);
%   s = dpmm_merge_samples(amdl, obs, H, ss, rthres);
%
%       merges multiple DPMM samples drawn from a MCMC chain (those
%       produced by the method make_output of dpmm class) into a single
%       sample by voting.
%
%       Input arguments:
%       - amdl:     the underlying atom model
%
%       - obs:      the associated observations
%
%       - H:        the struct of inheritance (empty if no inheritance)
%
%       - ss:       the sequence of samples
%
%       - rthres:   the threshold of ratio. If the count of an atom
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

n = amdl.get_num_samples(obs);

if ~isstruct(ss) && isvector(ss)
    error('dpmm_merge_samples:invalidarg', 'ss should be a struct vector.');
end

if ~isempty(H)
    if ~(isstruct(H) && isfield(H, 'tag') && strcmp(H.tag, 'dpmm_inherits'))
        error('dpmm_merge_samples:invalidarg', ...
            'H should be either empty or a dpmm_inherits struct.');
    end
    Kp = H.num;
else
    Kp = 0;
end

if nargin < 5
    rthres = 0;
else
    if ~(isfloat(rthres) && isscalar(rthres) && isreal(rthres) && ...
            rthres >= 0 && rthres < 1)
        error('dpmm_merge_samples:invalidarg', ...
            'rthres should be a non-negative scalar in [0, 1).');
    end
end


%% main

% get maximum atom id

max_aid = max([ss.max_atom_id]);

% voting

all_labels = vertcat(ss.labels);
labels = mode(all_labels, 1);
[atom_ids, counts, grps, z] = uniqvalues(labels, 'CGI');
K = numel(atom_ids);

% discard 

unlabeled = [];

if rthres > 0
    cthres = rthres * n;
    if any(counts < cthres)
        di = find(counts < cthres);
        if numel(di) == K
            error('dpmm_merge_samples:rterror', ...
                'No atom can be retained after thresholding.');
        end
        
        unlabeled = [grps{di}]; 
        
        atom_ids(di) = [];
        counts(di) = []; %#ok<NASGU>
        grps(di) = [];
        
        zmap = zeros(1, K);
        ri = setdiff(1:K, di);
        zmap(ri) = 1 : numel(ri);
        K = numel(ri);        
        z = zmap(z);
    end
end

% re-draw atoms 

atoms = cell(1, K);
if Kp == 0
    for k = 1 : K
        atoms{k} = amdl.posterior_atom(obs, grps{k});
    end
else
    [is_inherit, pids] = ismember(atom_ids, H.atom_ids);
    for k = 1 : K
        if is_inherit(k)
            atoms{k} = amdl.posterior_atom(obs, grps{k}, H.atoms{pids(k)});
        else
            atoms{k} = amdl.posterior_atom(obs, grps{k});
        end
    end
end

% relabel observations whose atoms were discarded

if ~isempty(unlabeled)
    
    assert(all(z(unlabeled) == 0));
    
    L = zeros(K, numel(unlabeled));    
    for k = 1 : K        
        llik = amdl.evaluate_loglik(atoms{k}, obs);
        L(k, :) = llik(unlabeled);
    end    
    [~, zu] = max(L, [], 1);
    
    z(unlabeled) = zu;
end

% make output

s.natoms = K;
s.max_atom_id = max_aid;
s.atom_ids = atom_ids;
s.atoms = atoms;
s.atom_counts = intcount(K, z);
s.iatoms = z;
s.labels = s.atom_ids(z);


