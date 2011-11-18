function modules = smi_modules()
% Get modules information of SMI toolbox
%
%   smi_modules();
%       prints information about modules of SMI toolbox.
%
%   modules = smi_modules();
%       get a list of module information.
%
%       modules is an array of structs, of which each element corresponds
%       to a particular module.
%
%       The struct has the following fields:
%       - name:         the name of the module.
%       - summary:      a brief description that summarizes the contents
%                       of the module.
%       - depends:      the dependent modules
%       - subpaths:     the path (relative to the root of smitoolbox)
%                       of each directory for this module.
%       - mex:          the list of mex targets   
%
%   Remarks
%   -------
%       This function reads all information from the configuration
%       file: smi_modules.cfg.
%

%   Created by Dahua Lin, on Nov 10, 2010
%

%% read file

cdir = fileparts(mfilename('fullpath'));
cfgfile = fullfile(cdir, 'smi_modules.cfg');

if ~exist(cfgfile, 'file')
    error('smi_modules:filenotfound', 'smi_modules.cfg is not found');
end

fid = fopen(cfgfile, 'r');
lines = textscan(fid, '%s', 'delimiter', '\n');
lines = lines{1};
fclose(fid);


%% read configuration

nm = 0;
cmodule = [];

n = length(lines);
i = 1;
while i <= n;
    
    % read current line
    line = strtrim(lines{i});    
    if isempty(line) || line(1) == '#'
        i = i + 1;
        continue;
    end
    
    % parse the current line
    
    if line(1) == '[' && line(end) == ']'  % module head
        
        % add current module to list
        if ~isempty(cmodule)
            nm = nm + 1;
            modules{nm, 1} = cmodule; %#ok<AGROW>
        end
        
        % start a new module
        name = strtrim(line(2:end-1));
        if isempty(name)
            parse_error(i, 'Invalid module head');
        end
        
        cmodule = struct( ...
            'name', name, ...
            'summary', [], ...
            'depends', [], ...
            'subpaths', [], ...
            'mex', []);
        
        i = i + 1;
        
    elseif line(1) == ':'  % summary
        
        [cmodule.summary, i] = read_summary(lines, i);
        
                
    elseif line(1) == '.'  % items starts
        
        [name, items, i] = read_items(lines, i);  
        
        switch (name)
            case {'depends', 'subpaths' }
                cmodule.(name) = items;
                
            case 'mex'
                nmex = numel(items);
                cmodule.mex = cell(nmex, 1);
                for j = 1 : nmex
                    cmodule.mex{j} = process_mex(items{j});
                end
                cmodule.mex = vertcat(cmodule.mex{:});
        end
        
    end
            
end


if ~isempty(cmodule)
    nm = nm + 1;
    modules{nm, 1} = cmodule;
end

modules = vertcat(modules{:});



%% sub functions for parsing

function [s, i] = read_summary(lines, i)

line = lines{i};
s{1} = strtrim(line(2:end));
n = 1;
i = i + 1;

while i < numel(lines)
    line = strtrim(lines{i});
    if ~isempty(line)
        n = n + 1;
        s{n} = [' ', line]; %#ok<AGROW>
        i = i + 1;
    else
        break;
    end
end

s = [s{:}];


function [name, items, i] = read_items(lines, i)

[name, remain] = strtok(lines{i}, '=');
name = strtrim(name(2:end));
remain = strtrim(remain(2:end));

if isempty(name) || isempty(remain)
    parse_error(i, 'Invalid line in starting an item.');
end

if remain(1) ~= '{'
    parse_error(i, 'Invalid line in starting an item.');
end
remain = remain(2:end);

if ~isempty(remain) && remain(end) == '}'
    remain = strtrim(remain(1:end-1)); 
    is_end = true;
else
    is_end = false;
end

rs = {remain};
while ~is_end
    i = i + 1;
    remain = strtrim(lines{i});    
    
    while (isempty(remain) || remain(1) == '#') && i < numel(lines)
        i = i + 1;
        remain = strtrim(lines{i});
    end
    
    if isempty(remain) 
        break;
    else                
        if remain(end) == '}'
            remain = strtrim(remain(1:end-1));
            is_end = true;
        else
            is_end = false;
        end
        rs = [rs; {remain}]; %#ok<AGROW>
    end
end
i = i + 1;


for k = 1 : numel(rs)
    r = rs{k};
    if ~isempty(r)
        toks = textscan(r, '%s', 'delimiter', ',');
        toks = toks{1};
    else
        toks = {};
    end
    rs{k} = toks;
end

items = vertcat(rs{:});

    
function r = process_mex(line)

[tname, remain] = strtok(line, ':');
tname = strtrim(tname);
remain = strtrim(remain(2:end));

if isempty(remain)
    error('smi_module:parseerror', 'Invalid line for mex target %s', tname);
end

srcs = textscan(remain, '%s').';

r.name = tname;
r.sources = srcs{1};


%% Auxiliary functions


function parse_error(i, msg)

error('smi_module:parseerror', ...
    'Parse error at line %d: %s', i, msg);

