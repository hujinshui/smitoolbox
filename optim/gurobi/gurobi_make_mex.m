function gurobi_mex(varargin)
% Perform mex compilation of Gurobi wrappers
%
%   gurobi_mex
%

% Created by Dahua Lin, on Jan 20, 2011
%


gdir = getenv('GUROBI_HOME');
if isempty(gdir)
    error('gurobi_mex:invalidarg', ...
        'The environment variable GUROBI_HOME is not set.');
end

incdir = fullfile(gdir, 'include');
libdir = fullfile(gdir, 'lib');

codedir = fullfile(fileparts(mfilename('fullpath')), 'private');

preargs = {'-outdir', codedir, ['-I' incdir], ['-L', libdir], '-lgurobi40', '-O'};
mex(preargs{:}, fullfile(codedir, 'gurobi_solve_mex.cpp'));




