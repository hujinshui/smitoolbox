function S = smi_allocate_struct(fnames, m, n)
% allocate a struct array with specified fieldnames 

nf = numel(fnames);
args = cell(2, nf);
args(1, :) = fnames;

S = repmat(struct(args{:}), m, n);
