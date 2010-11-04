function smi_build_mex(varargin)
% Re-build all mex files in the toolbox
%
%   smi_build_mex;
%   smi_build_mex name1 name2 ...
%

% Created by Dahua Lin, on Apr 7, 2010
%

%% main


% get environment 

boost_home = getenv('BOOST_HOME');
if isempty(boost_home)
    error('smi_build_mex:enverror', ...
        'Cannot find Boost C++ Library. Please add environment variable BOOST_HOME.');
end


% read in the build list

text = read_list('smi_mex_list.txt');

% parse the list

S = cellfun(@parse_entry, text);
n = length(S);

% filter the list
if ~isempty(varargin)
    [tf, si] = ismember(varargin, {S.name});
    if ~all(tf)
        error('smi_build_mex:invalidtarget', ...
            '%s is not a valid mex target.', ...
            varargin{find(~tf, 1)});
    end
    S = S(si);
    n = length(S);
end

% build 

bopts = {['-I' boost_home], '-O'};

rootdir = fileparts(fileparts(mfilename('fullpath')));

for i = 1 : n
    s = S(i);
    
    fprintf('Building %s ...\n', s.name);
        
    outdir = fullfile(rootdir, s.outdir);
    sources = cellfun(@(s) fullfile(rootdir, s), s.sources, 'UniformOutput', false);
    
    args = [bopts, {'-outdir', outdir}, sources];
    cmd = joinstr('mex', args{:});    
    fprintf('%s\n', cmd);
    fprintf('\n');
    
    mex(args{:});
    
    fprintf('\n');
end




%% parse function

function text = read_list(filename)

text = cell(10, 1);

fid = fopen(filename);
if fid < 0
    error('smi_build_mex:ioerror', 'Failed to open %s', filename);
end

i = 0;
while 1
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    
    line = strtrim(line);
    
    if ~isempty(line) && line(1) ~= '#'    
        i = i + 1;
        if i > numel(text)
            text{numel(text) * 2} = [];
        end
        text{i} = line;
    end
end

text = text(1:i);
fclose(fid);



function ps = parse_entry(line)

ts = regexp(line, '\s*:\s*', 'split');

if numel(ts) ~= 2 || isempty(ts{1}) || isempty(ts{2})
    error('smi_build_mex:parse_error', 'Invalid line %s', line);
end

name = ts{1};
srcs = regexp(ts{2}, '\s+', 'split');

ps.name = name;
ps.sources = srcs;
ps.outdir = fileparts(srcs{1});


function s = joinstr(varargin)

n = length(varargin);
terms = cell(1, 2 * n - 1);
terms(1:2:end) = varargin;
terms(2:2:end) = {' '};
s = [terms{:}];




