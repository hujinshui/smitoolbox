function S = dpmm_inherits(atom_ids, max_id, atoms, pricounts, q, T)
%DPMM_INHERITS Creates a struct that captures the inherited atoms for DPMM
%
%   S = DPMM_INHERITS(atom_ids, max_id, atoms, pricounts);
%   S = DPMM_INHERITS(atom_ids, max_id, atoms, pricounts, q);
%   S = DPMM_INHERITS(atom_ids, max_id, atoms, pricounts, q, T);
%
%       Creates a struct of priors based on inherited atoms for DPMM.
%
%       Input arguments:
%       - atom_ids:     the identifiers of inherited atoms
%       - atoms:        the cell array of inherited atoms
%       - pricounts:    the prior counts of the inherited atoms
%       - q:            the acceptance probabilities
%       - T:            the a function handle to transform atoms into
%                       either an atom or a distribution of atoms
%
%       If q is omitted, we set it to 1 by default.
%       If T is omitted, we set it to identity transform by default.
%
%       The output S is a struct with the following fields:
%       - tag:          a string: 'dpmm_inherits'
%       - num:          the number of atoms that can be inherited
%       - atom_ids:     the vector of inherited atom identifiers
%       - max_id:       the maximum value of previous id
%       - atoms:        the cell array of inherited (and transformed) atoms
%       - pricounts:    the vector of prior counts of the inherited atoms
%       - q:            the acceptance probability
%       

% Created by Dahua Lin, on Sep 25, 2011
%

%% verify input arguments

if ~(isnumeric(atom_ids) && isvector(atom_ids))
    error('dpmm_inherits:invalidarg', 'atom_ids should be a numeric vector.');
end

if ~(isnumeric(max_id) && isscalar(max_id) && max_id >= max(atom_ids))
    error('dpmm_inherits:invalidarg', 'max_id is invalid.');
end

if ~iscell(atoms)
    error('dpmm_inherits:invalidarg', 'atoms should be a cell array.');
end

if ~(isfloat(pricounts) && isvector(pricounts))
    error('dpmm_inherits:invalidarg', ...
        'pricounts should be a numeric vector.');
end

if ~isequal(size(atom_ids), size(atoms), size(pricounts))
    error('dpmm_inherits:invalidarg', ...
        'The sizes of atom_ids, atoms, and pricounts are inconsistent.');
end
na = numel(atoms);

if nargin >= 3
    if ~(isfloat(q) && isreal(q) && (isscalar(q) || numel(q) == na))
        error('dpmm_inherits:invalidarg', ...
            'q should be a real scalar or a vector of length == #atoms');
    end
    if isscalar(q)
        q = q * ones(1, na);
    else
        q = reshape(q, [1, na]);
    end
else
    q = zeros(1, na);
end

if nargin >= 4
    if ~isa(T, 'function_handle')
        error('dpmm_inherits:invalidarg', ...
            'T should be a function handle.');
    end
else
    T = [];
end

%% main

na = numel(atoms);

A = cell(1, na);
if isempty(T)
    for i = 1 : na
        A{i} = atoms{i};
    end
else
    for i = 1 : na
        A{i} = T(atoms{i});
    end
end

S.tag = 'dpmm_inherits';
S.num = na;
S.atom_ids = reshape(atom_ids, 1, na);
S.max_id = max_id;
S.atoms = A;
S.pricounts = pricounts;
S.q = q;

