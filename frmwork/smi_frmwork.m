classdef smi_frmwork < handle
    % The class to represent a statistical modeling and inference
    % framework
    %
    
    % Created by Dahua Lin, on Aug 24, 2011
    %
    
    
    %% properties
    
    % public properties 
    
    properties(GetAccess='public', SetAccess='private')
        
        num_vars = 0;       % the number of variables
        num_funcs = 0;      % the number of distinct sampling steps
        cycle_seq = [];     % the sequence of steps performed per cycle
                
        vars = [];          % the struct array comprised of variable specs        
        funcs = [];         % the struct array comprised of function specs
                                
        is_compiled = false;    % whether the frmwork has been compiled        
    end
    
    % private properties (implementation-specific)
    
    properties(GetAccess='public', SetAccess='private')
        
        var_map;        % the map: var name -> var id        
        func_map;       % the map: func name -> func id
        instructions;   % the compiled instructions
    end
    
    
    %% methods for frmwork building
    
    methods
        
        function add_var(frmwork, vname, vtype, vsize)
            % Add a named variable to the frmwork
            %
            %   frmwork.add_var(vname, vtype, vsize);
            %
            %       Adds a variable to the frmwork. Input arguments:
            %
            %       - vname:    the variable name
            %       - vtype:    the type of the variable
            %       - vsize:    the expected size of the variable
            %
            %   Note: the variable name must be a valid name of a
            %         matlab variable, i.e. isvarname(name) is true.
            %
            
            if frmwork.is_compiled
                error('smi_frmwork:invalidarg', ...
                    'Cannot add variable to a compiled frmwork.');
            end
            
            if ~isvarname(vname)
                error('smi_frmwork:invalidarg', ...
                    'The input vname is not a valid variable name.');
            end
            
            if ~(isvarname(vtype) && exist(vtype, 'class'))
                error('smi_frmwork:invalidarg', ...
                    'The input vtype is not a valid type name.');
            end
            
            vsize = smi_verify_vsize(vsize);
            if isempty(vsize)
                error('smi_frmwork:invalidarg', ...
                    'The input vsize is not a valid variable size spec.');
            end
                                        
            nv = frmwork.num_vars;
            new_var.name = vname;
            new_var.type = vtype;
            new_var.size = vsize;
            frmwork.vars = [frmwork.vars; new_var];
            frmwork.num_vars = nv+1;                   
        end
        
        
        function add_func(frmwork, name, func)
            % Add a named function to the frmwork and connect with variables
            %
            %   frmwork.add_func(name, func);
            %
            %       registers a named function to the frmwork that could 
            %       be used by sampling steps. The function should be
            %       an instance of a derived class of smi_func.
            %
            %       Note one can add function objects of the same
            %       class multiple times, provided that they are 
            %       with a different name, as each function object
            %       are uniquely identified with its name.
            %
            
            if frmwork.is_compiled
                error('smi_frmwork:invalidarg', ...
                    'Cannot add function to a compiled frmwork.');
            end
            
            if ~isvarname(name)
                error('smi_frmwork:invalidarg', ...
                    'The input name is not a valid step name.');
            end
            
            if ~isa(func, 'smi_func')
                error('smi_frmwork:invalidarg', ...
                    'The input func must be an smi_func object.');
            end            
                                                        
            nf = frmwork.num_funcs;            
            new_func.name = name;
            new_func.func = func;
            frmwork.funcs = [frmwork.funcs; new_func];
            frmwork.num_funcs = nf+1;            
        end
        
        
        function add_steps(frmwork, seq, nrepeat)
            % add one or multiple steps to each cycle
            %
            %   frmwork.add_steps(frmwork, step);
            %       appends one step to the end of each cycle
            %
            %   frmwork.add_steps(frmwork, step, nr);
            %       appends the input step, repeated nr times, to 
            %       the end of each cycle.
            %
            %   frmwork.add_steps(frmwork, seq);
            %       appends a sequence of steps, in form of a cell
            %       array of steps, to the end of each cycle.
            %
            %   frmwork.add_steps(frmwork, seq, nr);
            %       appends the input sequence of steps, repeated
            %       nr times as a whole, to the end of each cycle.
            %
            %   Each step is a char string in the following format:
            %
            %       <func_name> [slot-var-pairs]
            %
            %   Different slot-var-pairs are delimited by empty space,
            %   each pair is in form of <slot_name>:<var_name>.
            %   No space is allowed surrounding ':'.
            %
            
            if frmwork.is_compiled
                error('smi_frmwork:invalidarg', ...
                    'Cannot add step to a compiled frmwork.');
            end
            
            if ischar(seq) && ~isempty(seq) && ndims(seq) == 2 && size(seq,1) == 1
                seq = { seq };
            else
                if ~(iscellstr(seq) && isvector(seq))
                    error('smi_frmwork:invalidarg', ...
                        'seq must be a char string or a cell array of strings.');
                end    
                
                if size(seq, 2) > 1
                    seq = seq.';
                end
            end
                                 
            if nargin >= 3
                if ~(isnumeric(nrepeat) && isscalar(nrepeat) && ...
                        nrepeat == fix(nrepeat) && nrepeat >= 0)
                    error('smi_frmwork:invalidarg', ...
                        'nrepeat must be a non-negative integer scalar.');
                end
                
                if nrepeat == 0
                    return;
                end
            else
                nrepeat = 1;
            end
            
            ns = length(seq);
            ss = repmat( ...
                struct('func_name', [], 'slots', [], 'vars', []), ...
                ns, 1);
            
            for i = 1 : ns
                [ss(i).func_name, ss(i).slots, ss(i).vars] = ...
                    smi_frmwork.parse_step(seq{i});
            end
            
            if nrepeat > 1
                ss = repmat(ss, nrepeat, 1);
            end
            
            frmwork.cycle_seq = [frmwork.cycle_seq; ss];
        end           
    end
    
    
    methods(Static, Access='private')
        
        function [fname, slots, vs] = parse_step(statement)
            
            s = strtrim(statement);
            
            if isempty(s)
                error('smi_frmwork:parse_err', ...
                    'empty statement encountered .');
            end
            
            [fname, remain] = strtok(s);
            
            if ~isvarname(fname)
                error('smi_frmwork:parse_err', ...
                    'Invalid statement: %s', s);
            end
            
            ret = regexp(remain, '([\w\d]+):([\w\d]+)', 'tokens');
            np = numel(ret);
            
            slots = cell(np, 1);
            vs = cell(np, 1);
            
            for i = 1 : np
                slots{i} = ret{i}{1};
                vs{i} = ret{i}{2};
            end            
        end        
        
    end
        
    
    
    %% methods for displaying
    
    methods
        function disp(frmwork)
            % Display the basic information about the frmwork
            %
            %   disp(frmwork);
            %
            
            if frmwork.is_compiled
                fprintf('smi_frmwork [compiled]:\n');
            else
                fprintf('smi_frmwork [not compiled]:\n');
            end
            
            fprintf('--------------------------------\n');
            fprintf('\t# variables = %d\n', frmwork.num_vars);
            fprintf('\t# functions = %d\n', frmwork.num_funcs);
            fprintf('\t  cycle len = %d\n', length(frmwork.cycle_seq));
            fprintf('\n');
        end        
        
        
        function dump(frmwork, fid)
            % Display detailed information about a compiled frmwork
            %
            %   frmwork.dump();
            %       Prints the frmwork information to standard output.
            %
            %   frmwork.dump(fid);
            %       Prints the frmwork information to the file of
            %       the input file id.
            %                        
            
            if ~frmwork.is_compiled
                warning('smi_frmwork:noncompiledump', ...
                    'The smi_frmwork instance has not been compiled.');
                return;
            end
                                    
            % display variables
            
            if nargin < 2
                fid = 1;
            end
            
            fprintf(fid, 'Variable:\n');
            fprintf(fid, '------------\n');
            for i = 1 : frmwork.num_vars
                v = frmwork.vars(i);
                fprintf(fid, '  [%d] %s: %s %s\n', ...
                    i, v.name, v.type, smi_vsize2str(v.size));
            end
            fprintf(fid, '\n');
            
            % display functions                        
            
            fprintf(fid, 'Functions:\n');
            fprintf(fid, '------------\n');
            for i = 1 : frmwork.num_funcs
                f = frmwork.funcs(i);
                fo = f.func;
                fprintf(fid, '  [%d] %s: class %s\n', ...
                    i, f.name, class(fo));
                
                n_in = fo.num_input_slots;
                n_out = fo.num_output_slots;
                
                for j = 1 : n_in
                    slinfo = fo.get_slot_info('in', j);
                    fprintf(fid, '\t [in %d] %s: %s %s\n', ...
                        j, slinfo.name, slinfo.type, smi_vsize2str(slinfo.size));
                end
                
                for j = 1 : n_out
                    slinfo = fo.get_slot_info('out', j);
                    fprintf(fid, '\t [in %d] %s: %s %s\n', ...
                        j, slinfo.name, slinfo.type, smi_vsize2str(slinfo.size));
                end
            end
                                    
            % display cycle
            
            fprintf(fid, 'Cycle:\n');
            fprintf(fid, '------------\n');
            
            cseq = frmwork.cycle_seq;
            for i = 1 : length(cseq)
                ci = frmwork.instructions(i);
                f = frmwork.funcs(ci.f_id);
                
                fprintf(fid, '  (%d) %s:', i, f.name);
                
                for j = 1 : ci.n_in
                    if ci.v_in(j) > 0
                        fprintf(fid, ' %s', frmwork.vars(ci.v_in(j)).name);                        
                    end
                end
                
                fprintf(fid, ' ->');
                
                for j = 1 : ci.n_out
                    if ci.v_out(j) > 0
                        fprintf(fid, ' %s', frmwork.vars(ci.v_out(j)).name);
                    end
                end                
                fprintf(fid, '\n');
            end
            
            % blank line
            fprintf('\n');            
        end
        
    end
           
    
    %% methods for compilation
    
    methods
        
        function compile(frmwork)
            % Compile the frmwork
            %
            %   compile(frmwork) or frmwork.compile();
            %
            %       Compiles the frmwork such that it can be run more
            %       efficiently during sampling.
            %
            %       If the frmwork has alread been compiled, this 
            %       method is a no-action.
            %
            
            if frmwork.is_compiled
                return;
            end
           
            vmap = smi_frmwork.build_idmap({frmwork.vars.name});           
            fmap = smi_frmwork.build_idmap({frmwork.funcs.name});
            insts = smi_frmwork.compile_steps(frmwork.cycle_seq, ...
                frmwork.vars, vmap, frmwork.funcs, fmap); 
            
            frmwork.var_map = vmap;
            frmwork.func_map = fmap;
            frmwork.instructions = insts;

            frmwork.is_compiled = true;
        end                
    end
    
    
    methods(Static, Access='private')
        % private implementation of compilation
                
        function idmap = build_idmap(Vs, itemname)
            nv = numel(Vs);
            
            if numel(unique(Vs)) < nv
                error('smi_frmwork:comperr', ...
                    'Duplicated %s encountered.', itemname);
            end
            
            for i = 1 : nv
                idmap.(Vs{i}) = i;
            end
        end
        
        
        function insts = compile_steps(steps, V, vmap, F, fmap)
            % compile steps into instructions
            
            ns = numel(steps);
            insts = repmat( struct( ...
                'f_id', [], ...
                'n_in', [], ...
                'v_in', [], ...
                'n_out', [], ...
                'v_out', [], ...
                'out_flags', []), ...
                ns, 1);
            
            for i = 1 : ns
                
                fname = steps(i).func_name;
                ss = steps(i).slots;
                vs = steps(i).vars;
                
                if ~isfield(fmap, fname)
                    error('smi_frmwork:compile_err', ...
                        'Undeclared function name %s', fname);
                end
                f_id = fmap.(fname);
                fo = F(f_id).func;
                
                n_in = fo.num_input_slots;
                n_out = fo.num_output_slots;                                                
                v_in = zeros(1, n_in);
                v_out = zeros(1, n_out);
                                                    
                % set in/out map
                                               
                ni = numel(ss);
                for j = 1 : ni
                    [sdir, sid] = fo.get_slot_id(ss{j});
                    si = fo.get_slot_info(sdir, sid);
                    vname = vs{j};
                    if ~isfield(vmap, vname)
                        error('smi_frmwork:compile_err', ...
                            'Undeclared variable name %s', vname);
                    end
                    vid = vmap.(vname);
                    
                    if strcmpi(sdir, 'in')
                        v_in(sid) = vid;
                        if ~(smi_is_vcompatible(V(vid), si))
                            error('smi_frmwork:compile_err', ...
                                'Mismatch var %s => slot %s of function %s', ...
                                vs{j}, ss{j}, F(i).name);
                        end
                        
                    elseif strcmpi(sdir, 'out')
                        v_out(sid) = vid;
                        if ~(smi_is_vcompatible(si, V(vid)))
                            error('smi_frmwork:compile_err', ...
                                'Mismatch var %s <= slot %s of function %s', ...
                                vs{j}, ss{j}, F(i).name);
                        end
                        
                    end
                end                                
                
                % test whether slots are properly connected
                                
                if ~(fo.test_slots(v_in > 0, v_out > 0))
                    error('smi_frmwork:compile_err', ...
                        'The slots of function %s are not properly connected.', ...
                        F(i).name);
                end                                
                
                % test whether type matches
                                                                
                % store info
                insts(i).f_id = f_id;
                insts(i).n_in = n_in;              
                insts(i).v_in = v_in;                
                insts(i).n_out = n_out;
                insts(i).v_out = v_out;
                insts(i).out_flags = v_out > 0;
            end
            
        end
        
    end
       
end

