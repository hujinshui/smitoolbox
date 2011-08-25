classdef xgs_model < handle
    % The class to represent a probabilistic model for 
    % eXtensible Gibbs Sampling (XGS)
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
                                
        is_compiled = false;    % whether the model has been compiled        
    end
    
    % private properties (implementation-specific)
    
    properties(GetAccess='private', SetAccess='private')
        
        var_map;        % the map: var name -> var id        
        func_map;       % the map: func name -> func id
        instructions;   % the compiled instructions
    end
    
    
    %% methods for model building
    
    methods
        
        function add_var(model, vname, vtype, vsize)
            % Add a named variable to the model
            %
            %   model.add_var(vname, vtype, vsize);
            %
            %       Adds a variable to the model. Input arguments:
            %
            %       - vname:    the variable name
            %       - vtype:    the type of the variable
            %       - vsize:    the expected size of the variable
            %
            %   Note: the variable name must be a valid name of a
            %         matlab variable, i.e. isvarname(name) is true.
            %
            
            if model.is_compiled
                error('xgs_model:invalidarg', ...
                    'Cannot add variable to a compiled model.');
            end
            
            if ~isvarname(vname)
                error('xgs_model:invalidarg', ...
                    'The input vname is not a valid variable name.');
            end
            
            if ~(isvarname(vtype) && exist(vtype, 'class'))
                error('xgs_model:invalidarg', ...
                    'The input vtype is not a valid type name.');
            end
            
            vsize = xgs_verify_vsize(vsize);
            if isempty(vsize)
                error('xgs_model:invalidarg', ...
                    'The input vsize is not a valid variable size spec.');
            end
                                        
            nv = model.num_vars;
            new_var.name = vname;
            new_var.type = vtype;
            new_var.size = vsize;
            model.vars = [model.vars; new_var];
            model.num_vars = nv+1;                   
        end
        
        
        function add_func(model, name, func)
            % Add a named function to the model and connect with variables
            %
            %   model.add_func(name, func);
            %
            %       registers a named function to the model that could 
            %       be used by sampling steps. The function should be
            %       an instance of a derived class of xgs_func.
            %
            %       Note one can add function objects of the same
            %       class multiple times, provided that they are 
            %       with a different name, as each function object
            %       are uniquely identified with its name.
            %
            
            if model.is_compiled
                error('xgs_model:invalidarg', ...
                    'Cannot add function to a compiled model.');
            end
            
            if ~isvarname(name)
                error('xgs_model:invalidarg', ...
                    'The input name is not a valid step name.');
            end
            
            if ~isa(func, 'xgs_func')
                error('xgs_model:invalidarg', ...
                    'The input func must be an xgs_func object.');
            end            
                                                        
            nf = model.num_funcs;            
            new_func.name = name;
            new_func.func = func;
            model.funcs = [model.funcs; new_func];
            model.num_funcs = nf+1;            
        end
        
        
        function add_steps(model, seq, nrepeat)
            % add one or multiple steps to each cycle
            %
            %   model.add_steps(model, step);
            %       appends one step to the end of each cycle
            %
            %   model.add_steps(model, step, nr);
            %       appends the input step, repeated nr times, to 
            %       the end of each cycle.
            %
            %   model.add_steps(model, seq);
            %       appends a sequence of steps, in form of a cell
            %       array of steps, to the end of each cycle.
            %
            %   model.add_steps(model, seq, nr);
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
            
            if model.is_compiled
                error('xgs_model:invalidarg', ...
                    'Cannot add step to a compiled model.');
            end
            
            if ischar(seq) && ~isempty(seq) && ndims(seq) == 2 && size(seq,1) == 1
                seq = { seq };
            else
                if ~(iscellstr(seq) && isvector(seq))
                    error('xgs_model:invalidarg', ...
                        'seq must be a char string or a cell array of strings.');
                end    
                
                if size(seq, 2) > 1
                    seq = seq.';
                end
            end
                                 
            if nargin >= 3
                if ~(isnumeric(nrepeat) && isscalar(nrepeat) && ...
                        nrepeat == fix(nrepeat) && nrepeat >= 0)
                    error('xgs_model:invalidarg', ...
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
                    xgs_model.parse_step(seq{i});
            end
            
            if nrepeat > 1
                ss = repmat(ss, nrepeat, 1);
            end
            
            model.cycle_seq = [model.cycle_seq; ss];
        end           
    end
    
    
    methods(Static, Access='private')
        
        function [fname, slots, vs] = parse_step(statement)
            
            s = strtrim(statement);
            
            if isempty(s)
                error('xgs_model:parse_err', ...
                    'empty statement encountered .');
            end
            
            [fname, remain] = strtok(s);
            
            if ~isvarname(fname)
                error('xgs_model:parse_err', ...
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
        function disp(model)
            % Display the basic information about the model
            %
            %   disp(model);
            %
            
            if model.is_compiled
                fprintf('xgs_model [compiled]:\n');
            else
                fprintf('xgs_model [not compiled]:\n');
            end
            
            fprintf('------------------------\n');
            fprintf('\t# variables = %d\n', model.num_vars);
            fprintf('\t# functions = %d\n', model.num_funcs);
            fprintf('\t  cycle len = %d\n', length(model.cycle_seq));
            fprintf('\n');
        end        
        
        
        function dump(model, fid)
            % Display detailed information about a compiled model
            %
            %   model.dump();
            %       Prints the model information to standard output.
            %
            %   model.dump(fid);
            %       Prints the model information to the file of
            %       the input file id.
            %                        
            
            if ~model.is_compiled
                warning('xgs_model:noncompiledump', ...
                    'The xgs_model has not been compiled.');
                return;
            end
                                    
            % display variables
            
            if nargin < 2
                fid = 1;
            end
            
            fprintf(fid, 'Variable:\n');
            fprintf(fid, '------------\n');
            for i = 1 : model.num_vars
                v = model.vars(i);
                fprintf(fid, '  [%d] %s: %s %s\n', ...
                    i, v.name, v.type, xgs_vsize2str(v.size));
            end
            fprintf(fid, '\n');
            
            % display functions                        
            
            fprintf(fid, 'Functions:\n');
            fprintf(fid, '------------\n');
            for i = 1 : model.num_funcs
                f = model.funcs(i);
                fo = f.func;
                fprintf(fid, '  [%d] %s: class %s\n', ...
                    i, f.name, class(fo));
                
                n_in = fo.num_input_slots;
                n_out = fo.num_output_slots;
                
                for j = 1 : n_in
                    slinfo = fo.get_slot_info('in', j);
                    fprintf(fid, '\t [in %d] %s: %s %s\n', ...
                        j, slinfo.name, slinfo.type, xgs_vsize2str(slinfo.size));
                end
                
                for j = 1 : n_out
                    slinfo = fo.get_slot_info('out', j);
                    fprintf(fid, '\t [in %d] %s: %s %s\n', ...
                        j, slinfo.name, slinfo.type, xgs_vsize2str(slinfo.size));
                end
            end
                                    
            % display cycle
            
            fprintf(fid, 'Cycle:\n');
            fprintf(fid, '------------\n');
            
            cseq = model.cycle_seq;
            for i = 1 : length(cseq)
                ci = model.instructions(i);
                f = model.funcs(ci.f_id);
                
                fprintf(fid, '  (%d) %s:', i, f.name);
                
                for j = 1 : ci.n_in
                    if ci.v_in(j) > 0
                        fprintf(fid, ' %s', model.vars(ci.v_in(j)).name);                        
                    end
                end
                
                fprintf(fid, ' ->');
                
                for j = 1 : ci.n_out
                    if ci.v_out(j) > 0
                        fprintf(fid, ' %s', model.vars(ci.v_out(j)).name);
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
        
        function compile(model)
            % Compile the model
            %
            %   compile(model) or model.compile();
            %
            %       Compiles the model such that it can be run more
            %       efficiently during sampling.
            %
            %       If the model has alread been compiled, this 
            %       method is a no-action.
            %
            
            if model.is_compiled
                return;
            end
           
            vmap = xgs_model.build_idmap({model.vars.name});           
            fmap = xgs_model.build_idmap({model.funcs.name});
            insts = xgs_model.compile_steps(model.cycle_seq, ...
                model.vars, vmap, model.funcs, fmap); 
            
            model.var_map = vmap;
            model.func_map = fmap;
            model.instructions = insts;

            model.is_compiled = true;
        end                
    end
    
    
    methods(Static, Access='private')
        % private implementation of compilation
                
        function idmap = build_idmap(Vs, itemname)
            nv = numel(Vs);
            
            if numel(unique(Vs)) < nv
                error('xgs_model:comperr', ...
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
                    error('xgs_model:compile_err', ...
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
                        error('xgs_model:compile_err', ...
                            'Undeclared variable name %s', vname);
                    end
                    vid = vmap.(vname);
                    
                    if strcmpi(sdir, 'in')
                        v_in(sid) = vid;
                        if ~(xgs_is_vcompatible(V(vid), si))
                            error('xgs_model:compile_err', ...
                                'Mismatch var %s => slot %s of function %s', ...
                                vs{j}, ss{j}, F(i).name);
                        end
                        
                    elseif strcmpi(sdir, 'out')
                        v_out(sid) = vid;
                        if ~(xgs_is_vcompatible(si, V(vid)))
                            error('xgs_model:compile_err', ...
                                'Mismatch var %s <= slot %s of function %s', ...
                                vs{j}, ss{j}, F(i).name);
                        end
                        
                    end
                end                                
                
                % test whether slots are properly connected
                                
                if ~(fo.test_slots(v_in > 0, v_out > 0))
                    error('xgs_model:compile_err', ...
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
    
    
    %% methods for run
    
    methods
        
        function S = run(model, opts, initcfg)
            % Run a simulation procedure from an initial sample
            %
            %   S = model.run(opts, initsample);
            %
            %       Runs the simulation, starting from the provided
            %       initial sample.
            %
            %       Arguments:
            %       - opts:     the options that controls the sampling
            %                   procedure. 
            %                   It can be obtained by calling xgs_options.
            %
            %       - initcfg:  the initial configuration to begin with.
            %
            %       - S:        A struct array comprised of obtained
            %                   samples. S(i) is the i-th sample. 
            %                   The fields are the variable names.
            %            
                                                
            opts = xgs_options(opts);
            
            if ~(isstruct(initcfg) && isscalar(initcfg))
                error('xgs_model:invalidarg', ...
                    'initcfg should be a struct scalar.');
            end
            repr = preprocess_sample(model, initcfg);
            
            S = do_run(model, opts, repr, 0);            
        end
        
        function R = par_run(model, opts, initcfgs)
            % Run parallel simulations
            %
            %   S = model.run(opts, initsample);
            %
            %       Runs the simulation, starting from multiple
            %       initial samples in parallel.
            %
            %       Arguments:
            %       - opts:     the options that controls the sampling
            %                   procedure. 
            %                   It can be obtained by calling xgs_options.
            %
            %       - initcfgs:  the initial configurations to begin with.
            %
            %       - R:        A cell array containing resultant samples.
            %                   In particular, R{i} is a struct array
            %                   comprised of the samples obtained through
            %                   the simulation starting from the i-th
            %                   initial sample.
            %   
            
            opts = xgs_options(opts);
                        
            if opts.display > 2
                warning('xgs_model:displayadjust', ...
                    'opts.display is adjusted to 2 for par_run.');
                opts.display = 2;
            end
            
            if ~isstruct(initcfgs)
                error('xgs_model:invalidarg', ...
                    'initcfgs should be a struct array.');
            end
            nc = numel(initcfgs); 
            
            C = cell(size(initcfgs));            
            parfor i = 1 : nc
                C{i} = preprocess_sample(model, initcfgs(i));
            end
            
            R = cell(size(initcfgs));
            parfor i = 1 : nc
                R{i} = do_run(model, opts, C{i}, i);
            end
        end
                
    end
    
    
    methods(Access='private')
        
        function repr = preprocess_sample(model, cfg)
            % form internal representation and verify
            
            nv = model.num_vars;
            vs = model.vars;
            
            % form internal representation
            repr = cell(1, nv);
            for i = 1 : nv
                vn = vs(i).name;
                if isfield(cfg, vn)
                    repr{i} = cfg.(vn);
                end
            end                       
        end
        
        
        function S = do_run(model, opts, repr, ithread)
            % run with a single config
                 
            % extract fields
            
            cseq = model.cycle_seq;
            clen = numel(cseq);
            
            fs = model.funcs;
            insts = model.instructions;
            vmap = model.var_map;
            
            T0 = opts.burnin;
            N = opts.nsamples;
            nc = opts.cps;
            displevel = opts.display;
            
            if displevel >= 1
                if ithread > 0
                    dprefix = ['[[thread ' int2str(ithread), ']] '];
                else
                    dprefix = '';
                end
            end
            
            % get ids of variable to be output to samples
            if isempty(opts.output_vars)
                ovars = {model.vars.name};
            else
                ovars = opts.output_vars;
            end
            ovar_ids = zeros(size(ovars));
            for i = 1 : numel(ovars)
                ovar_ids(i) = vmap.(ovars{i});
            end
            
            % initialize result container
            S = cell(1, N);
            
            % main loop
            Tstart = -T0;
            Tend = (nc - 1) * N;
            
            if displevel >= 2
                disp([dprefix, 'simulation started']);
            end
            
            c = 0; % sample counter            
            for t = Tstart : Tend
                
                if t == 0 && displevel >= 2
                    disp([dprefix, 'burn-in finished']);
                end
                
                if displevel >= 3                    
                    disp([dprefix, sprintf('  cycle %d ...', t)]);
                end
                                              
                % do update cycle
                for i = 1 : clen 
                                        
                    % pick the current step
                    st = insts(i);
                    f = fs(st.f_id).func;                  
                    
                    if displevel >= 4
                        disp([dprefix, sprintf('    step %d: %s', ...
                            i, fs(st.f_id).name)]);
                    end
                    
                    % form inputs                                        
                    vin = xgs_from_repr(repr, st.v_in);
                    
                    % do update computation
                    vout = cell(1, st.n_out);
                    [vout{:}] = f.evaluate(st.out_flags, vin{:});
                    
                    % put the results back to repr
                    for j = 1 : st.n_out
                        if st.v_out(j) > 0
                            repr{st.v_out(j)} = vout{j};
                        end
                    end
                end       
                                                
                % collect sample                
                if t >= 0 && mod(t, nc) == 0
                    c = c + 1;                    
                    for i = 1 : numel(ovars)
                        sp.(ovars{i}) = repr{ovar_ids(i)};
                    end
                    S{c} = sp;
                end
            end
            
            % convert results to struct form
            S = vertcat(S{:});  
            
            if displevel >= 1
                disp([dprefix, 'simulation done']);
            end
            
        end
        
    end

    
end

