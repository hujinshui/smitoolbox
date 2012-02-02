function info = smi_test(module, suitename)
% Performs unit testing for SMI modules
%
%   info = smi_test(module, suitename);
%   info = smi_test(module, suitenames);
%       Runs a test suite or a collection of suites under a module.
%
%       Here, module and suitename are both strings, and suitenames
%       is a cell array of strings.
%
%       In the output, info is a struct (or struct array, if for
%       a collection of multiple suites), which has the following field:
%
%       - 'module':         the module name
%       - 'name':           the suite name
%       - 'classname':      the class name
%       - 'ntotal':         the total number of cases
%       - 'nfailed':        the number of failed cases
%       - 'failed_cases':   the cell array of names of failed cases
%       
%   names = smi_test(module, 'get.suitenames');
%       Gets all test suite names under a module.   
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 3, 2011
%


%% verify input arguments

if ~isvarname(module)
    error('smi_test:invalidarg', 'The module name is not valid.');
end

snames_all = get_suitenames(module);

if nargin < 2
    use_all = 1;    
    snames = snames_all;    
    get_attr = [];
else
    use_all = 0;
    if ischar(suitename)
        if strcmpi(suitename, 'get.suitenames');
            get_attr = 'suitenames';
            snames = [];
        else
            snames = {suitename};
            get_attr = [];
        end
    elseif iscellstr(suitename)
        snames = suitename;
        get_attr = [];
    else
        error('smi_test:invalidarg', ...
            'The suitename should be either a string or a cell array of strings.');
    end
end

nsuites = numel(snames);

%% main

if ~isempty(get_attr)
    info = snames_all;
    return;
end

if ~use_all
    for i = 1 : nsuites
        sname = snames{i};
        if ~ismember(sname, snames_all)
            error('smi_test:invalidarg', ...
                'There is no test suite named %s for module %s.', ...
                sname, module);
        end
    end
end

info.module = module;
info.name = [];
info.classname = [];
info.ntotal = 0;
info.nfailed = 0;
info.failed_cases = {};

if nsuites > 1
    info = repmat(info, [nsuites, 1]);
end

fid = 1;

for i = 1 : nsuites

    sname = snames{i};
    
    fprintf(fid, 'Test suite [%d / %d]: %s\n', i, nsuites, snames{i});
    fprintf(fid, '====================================\n');
    
    clsname = ['tsuite_' sname];
    
    tobj = feval(clsname);
    
    Ms = methods(tobj);
    cases = cell(numel(Ms), 1);
    nc = 0;
    
    for j = 1 : numel(Ms)
        cmethod = Ms{j};
        if length(cmethod) >= 6 && strcmp(cmethod(1:5), 'test_')
            nc = nc + 1;
            cases{nc} = cmethod;
        end
    end
    
    fcases = {};
    
    for j = 1 : nc
        cmethod = cases{j};
        fprintf(fid, '  running %s ...\n', cmethod);
        try
            feval(cmethod, tobj);
        catch err
            fcases = [fcases; {cmethod}]; 
            print_error(fid, cmethod, err);
        end
    end
    fprintf('  --------------------------------\n');
    if isempty(fcases)
        fprintf('  All %d cases passed\n', nc);
    else
        fprintf('  In all %d cases, %d failed\n', nc, numel(fcases));
    end
        
    fprintf(fid, '\n');
    
    % store information
    
    info(i).name = sname; %#ok<*AGROW>
    info(i).classname = clsname;
    info(i).ntotal = nc;
    info(i).nfailed = numel(fcases);
    info(i).failed_cases = fcases;
    
end



%% Auxiliary functions

function sall = get_suitenames(module)
% retrieve suite names

smiroot = fileparts(fileparts(mfilename('fullpath')));
moduleroot = fullfile(smiroot, module);

if exist(moduleroot, 'dir')
    
    testdir = fullfile(moduleroot, 'tests');
    
    sall = {};
    if exist(testdir, 'dir')
        fns = dir(fullfile(testdir, 'tsuite_*.m'));        
        if ~isempty(fns)
            ns = numel(fns);
            sall = cell(ns, 1);
            for i = 1 : numel(fns)
                [pstr, cname] = fileparts(fns(i).name);  %#ok<ASGLU>
                if ~exist(cname, 'class')
                    error('smi_test:rterror', ...
                        '%s is not a class or is not in search-path.', cname);
                end
                sall{i} = cname(8:end);
            end
        end
    end
    
    if isempty(sall)        
        warning('smi_test:notests', ...
            'There are no testes for module %s', module);
    end
        
else
    error('smi_test:nonexist', ...
        'Cannot find the module named %s', module);
end


function print_error(fid, cmethod, err)
% print an error

fprintf(fid, '    %s failed: %s\n', cmethod, err.message);
fprintf(fid, '    Stack: \n');
m = numel(err.stack) - 1;   % ignore the smi_test function
for i = 1 : m
    s = err.stack(i);
    [pstr, fname, fext] = fileparts(s.file);  %#ok<ASGLU>
    fprintf(fid, '      %s (line %d)\n', [fname fext], s.line);
end


