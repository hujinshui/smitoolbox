function solver = svm_default_solver()
% Returns the default QP solver for SVM training
%
%   solver = svm_default_solver();
%

optim_ver = ver('optim');
optim_ver = str2double(optim_ver.Version);

if optim_ver >= 6
    alg = 'interior-point-convex';
else
    alg = 'interior-point';
end
opts = optimset('Algorithm', alg, 'Display', 'off');

solver = @(pb) mstd_solve(pb, opts);
