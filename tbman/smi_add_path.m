function smi_add_path()
% The function to add smitoolbox to MATLAB search paths
%
%   smi_add_path;
%
%   Note that this function does not save path permantly. To do so,
%   one can call savepath.
%

%% main

mdls = smi_modules;

rdir = fileparts(fileparts(mfilename('fullpath')));

subpaths = ['tbman'; vertcat(mdls.subpaths)];
fpaths = cellfun(@(p) fullfile(rdir, p), subpaths, 'UniformOutput', false);

addpath(fpaths{:});

