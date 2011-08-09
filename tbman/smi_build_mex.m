function smi_build_mex(varargin)
% Re-build all mex files in the toolbox
%
%   smi_build_mex;
%   smi_build_mex name1 name2 ...
%


%% main

mdls = smi_modules;

bopts = {'-O'};

% get environment 

bcslib_home = getenv('BCSLIB_HOME');
if isempty(bcslib_home)
    error('smi_build_mex:enverror', ...
        'Cannot find BCS Lib. Please add environment variable BCSLIB_HOME.');
end
bopts = [bopts, {['-I' bcslib_home]}];


% require_boost = ismember('graph', {mdls.name});
% 
% if require_boost
%     boost_home = getenv('BOOST_HOME');
%     if isempty(boost_home)
%         error('smi_build_mex:enverror', ...
%             'Cannot find Boost C++ Library. Please add environment variable BOOST_HOME.');
%     end
%     
%     bopts = [bopts, {['-I' boost_home]}];
% end


% get the list

S = vertcat(mdls.mex);

% filter the list
if ~isempty(varargin)
    [tf, si] = ismember(varargin, {S.name});
    if ~all(tf)
        error('smi_build_mex:invalidtarget', ...
            '%s is not a valid mex target.', ...
            varargin{find(~tf, 1)});
    end
    S = S(si);
end
n = length(S);


% build 

rootdir = fileparts(fileparts(mfilename('fullpath')));

for i = 1 : n
    s = S(i);
    
    fprintf('Building %s ...\n', s.name);
            
    outdir = fullfile(rootdir, fileparts(s.sources{1}));
    sources = cellfun(@(s) fullfile(rootdir, s), s.sources, 'UniformOutput', false);
    
    args = [bopts, {'-outdir', outdir}, sources];
    cmd = joinstr('mex', args{:});    
    fprintf('%s\n', cmd);
    fprintf('\n');
    
    mex(args{:});
    
    fprintf('\n');
end


%% parse function


function s = joinstr(varargin)

n = length(varargin);
terms = cell(1, 2 * n - 1);
terms(1:2:end) = varargin;
terms(2:2:end) = {' '};
s = [terms{:}];




