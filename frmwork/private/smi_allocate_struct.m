function S = smi_allocate_struct(fnames, m, n)
% allocate a struct array with specified fieldnames 

if nargin < 2
    m = 1;
end

if nargin < 3
    n = 1;
end

nf = numel(fnames);
args = cell(2, nf);
args(1, :) = fnames;

S = struct(args{:});
if m ~= 1 || n ~= 1
    S = repmat(S, m, n);
end
