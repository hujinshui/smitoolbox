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
        num_blocks = 0;     % the number of blocks
        
        vars = cell(16, 1);     % the list of added variables
        funcs = cell(16, 1);    % the list of added functions
        blocks = cell(16, 1);   % the list of added blocks 
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
            
            if ~isvarname(vname)
                error('smi_frmwork:invalidarg', ...
                    'The input vname is not a valid variable name.');
            end
            
            if ~smi_is_valid_typename(vtype)
                error('smi_frmwork:invalidarg', ...
                    'The input vtype is not a valid type name.');
            end
            
            vsize = smi_verify_vsize(vsize);
            if isempty(vsize)
                error('smi_frmwork:invalidarg', ...
                    'The input vsize is not a valid variable size spec.');
            end
            
            new_var.name = vname;
            new_var.type = vtype;
            new_var.size = vsize;
            
            nv = frmwork.num_vars;
            capa = numel(frmwork.vars);
            if capa == nv
                frmwork.vars = [frmwork.vars; cell(capa, 1)];
            end
            nv = nv + 1;
            frmwork.vars{nv} = new_var;
            frmwork.num_vars = nv;
            
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

            if ~isvarname(name)
                error('smi_frmwork:invalidarg', ...
                    'The input name is not a valid step name.');
            end
            
            if ~isa(func, 'smi_func')
                error('smi_frmwork:invalidarg', ...
                    'The input func must be an smi_func object.');
            end            
                                                                  
            new_func.name = name;
            new_func.obj = func;
            
            nf = frmwork.num_funcs;
            capa = numel(frmwork.funcs);
            if capa == nf
                frmwork.funcs = [frmwork.funcs; cell(capa, 1)];
            end
            nf = nf + 1;            
            frmwork.funcs{nf} = new_func;
            frmwork.num_funcs = nf;
        end
        
        
        function add_step(frmwork, func_name, inputs, outputs, nr, host_blk)
            % Add a step to a block in the framework
            %
            %   frmwork.add_step(func_name, inputs, outputs);
            %   frmwork.add_step(func_name, inputs, outputs, nr);
            %   frmwork.add_step(func_name, inputs, outputs, nr, blk_name)
            %
            %       adds a step to a block in the framework. By default,
            %       the step is added to the default block
            %                   
            %       Input arguments:
            %       - func_name:    the name of the function to be called 
            %                       by this step
            %
            %       - inputs/outputs:   an m x 2 cell array, with the 
            %                           first column for slot names, and
            %                           the second column for variable
            %                           names
            %
            %       - nr:           the number of times to repeat this
            %                       action. (by default, nr = 1).
            %
            %       - host_blk:     the name of the block that the step
            %                       is added to.
            %
            
            % verify input arguments
            
            if ~isvarname(func_name)
                error('smi_frmwork:invalidarg', ...
                    'The input func_name is invalid.');
            end
            
            if ~isempty(inputs)
                if ~(iscell(inputs) && ndims(inputs) == 2 && size(inputs, 2) == 2)
                    error('smi_frmwork:invalidarg', ...
                        'inputs should be a cell array with two columns.');
                end
                if ~all(cellfun(@isvarname, inputs(:)))
                    error('smi_frmwork:invalidarg', ...
                        'Some names in inputs are not valid.');
                end
            end
            
            if ~isempty(outputs)
                if ~(iscell(outputs) && ndims(outputs) == 2 && size(outputs, 2) == 2)
                    error('smi_frmwork:invalidarg', ...
                        'outputs should be a cell array with two columns.');
                end
                if ~all(cellfun(@isvarname, outputs(:)))
                    error('smi_frmwork:invalidarg', ...
                        'Some names in outputs are not valid.');
                end
            end                        
            
            if nargin < 5
                nr = 1;
            else
                if ~(isnumeric(nr) && isscalar(nr) && nr == fix(nr) && nr >= 0)
                    error('smi_frmwork:invalidarg', ...
                        'nr should be a non-negative integer scalar.');
                end
            end
            
            if nargin < 6
                host_blk = 'default';
            else
                if ~isvarname(host_blk)
                    error('smi_frmwork:invalidarg', 'host_blk is invalid.');
                end
            end
                        
            % create and add new step
            
            new_step.type = 'f';
            new_step.func_name = func_name;
            new_step.nrepeats = nr;
            
            if ~isempty(inputs)
                new_step.input_slots = inputs(:, 1);
                new_step.input_vars = inputs(:, 2);
            else
                new_step.input_slots = {};
                new_step.input_vars = {};
            end
            
            if ~isempty(outputs)
                new_step.output_slots = outputs(:, 1);
                new_step.output_vars = outputs(:, 2);
            else
                new_step.output_slots = {};
                new_step.output_vars = {};
            end
            
            do_add_step(frmwork, new_step, host_blk);
                       
        end
        
        
        function add_block_step(frmwork, blk_name, nr, host_blk)
            % Add a block as a step in another block
            %
            %   frmwork.add_block_step(frmwork, blk_name, nr, host_blk);
            %
            %       adds the specified block as a step in the host block.
            %
            %       Input arguments:
            %       - blk_name:     the name of the block to be added
            %       - nr:           the number of times to repeat
            %       - host_blk:     the name of the host block
            %
            
            % verify input
            
            if ~isvarname(blk_name)
                error('smi_frmwork:invalidarg', 'blk_name is invalid.');
            end
            
            if nargin < 3
                nr = 1;
            else
                if ~(isnumeric(nr) && isscalar(nr) && nr == fix(nr) && nr >= 0)
                    error('smi_frmwork:invalidarg', ...
                        'nr should be a non-negative integer scalar.');
                end
            end
            
            if nargin < 4
                host_blk = 'default';
            else
                if ~isvarname(host_blk)
                    error('smi_frmwork:invalidarg', 'host_blk is invalid.');
                end
            end            
            
            new_step.type = 'b';
            new_step.block_name = blk_name;
            new_step.nrepeats = nr;
            
            do_add_step(frmwork, new_step, host_blk);            
        end
        
                
    end
    
    
    methods(Access='private')
        
        function do_add_step(frmwork, new_step, blk_name)
            
            % search or create the block            
            nb = frmwork.num_blocks;
            B = frmwork.blocks;
            iblk = 0;
            for i = 1 : nb
                if strcmp(B{i}.name, blk_name)
                    iblk = i;
                    break;
                end
            end
            
            if iblk == 0 % not found, and create new block
                new_blk.name = blk_name;
                new_blk.num_steps = 0;
                new_blk.steps = cell(16, 1);
                
                capa = numel(frmwork.blocks);
                if capa == nb
                    frmwork.blocks = [frmwork.blocks; cell(capa, 1)];
                end
                iblk = nb + 1;
                frmwork.blocks{iblk} = new_blk;
                frmwork.num_blocks = iblk;
            end
            
            % add step to the block
            
            ns = frmwork.blocks{iblk}.num_steps;
            capa = numel(frmwork.blocks{iblk}.steps);
            
            if capa == ns
                frmwork.blocks{iblk}.steps = ... 
                    [frmwork.blocks{iblk}.steps; cell(capa, 1)];
            end
            
            ns = ns + 1;
            frmwork.blocks{iblk}.steps{ns} = new_step; 
            frmwork.blocks{iblk}.num_steps = ns;
        end
        
    end
    
    
    %% methods for displaying
    
    methods
        
        function disp(frmwork)
            % Display the framework
            %
            %   disp(frmwork);
            %
            
            fprintf('SMI Framework\n');
            fprintf('------------------------\n');
            
            fprintf('  # Variables = %d\n', frmwork.num_vars);
            fprintf('  # Functions = %d\n', frmwork.num_funcs);
            fprintf('  # Blocks = %d\n', frmwork.num_blocks);
            fprintf('\n');
        end
        
        
        function dump(frmwork, fid)
            % Display the framework in detail
            %
            %   dump(frmwork);
            %       Displays detailed information of the framework to
            %       console
            %
            %   dump(frmwork, fid);
            %       Prints the detailed information of the framework 
            %       to a file (with given file id).
            %
            
            if nargin < 2
                fid = 1;
            end
            
            fprintf(fid, 'SMI Framework (Details)\n');
            fprintf(fid, '===============================\n');
            
            fprintf(fid, 'Variables\n');
            fprintf(fid, '----------------\n');
            
            nv = frmwork.num_vars;
            for i = 1 : nv
                v = frmwork.vars{i};
                smi_frmwork.dump_var(fid, 1, v);
            end
            fprintf(fid, '\n');
            
            fprintf(fid, 'Functions\n');
            fprintf(fid, '----------------\n');
            
            nf = frmwork.num_funcs;
            for i = 1 : nf
                f = frmwork.funcs{i};
                fo = f.obj;
                fprintf(fid, '    %s: class %s\n', f.name, class(fo)); 
            end
            fprintf(fid, '\n');
            
            fprintf(fid, 'Blocks\n');
            fprintf(fid, '----------------\n');
            
            nb = frmwork.num_blocks;
            for i = 1 : nb
                b = frmwork.blocks{i};
                ns = b.num_steps;
                fprintf(fid, 'Block %s (#steps = %d)\n', b.name, ns);
                                
                for j = 1 : ns
                    s = b.steps{j};
                    if s.type == 'b'
                        fprintf(fid, '    [%d] call block %s [#=%d]\n', ...
                            j, s.block_name, s.nrepeats);
                    elseif s.type == 'f'
                        fprintf(fid, '    [%d] call func %s [#=%d]: (%s) => (%s)\n', ...
                            j, s.func_name, s.nrepeats, ...
                            smi_frmwork.conns2str(s.input_slots, s.input_vars), ...
                            smi_frmwork.conns2str(s.output_slots, s.output_vars));
                    end
                end
            end
            fprintf(fid, '\n');
            
        end
        
    end
    
    
    
    methods(Static, Access='private')
        
        function dump_var(fid, tabs, v)            
            prefix = blanks(4 * tabs);
            fprintf(fid, [prefix, '%s: %s %s\n'], ...
                v.name, v.type, smi_vsize2str(v.size)); 
        end
        
        
        function s = conns2str(slots, vars)
            n = numel(slots);
            if n > 0
                terms = cell(1, n);
                for i = 1 : n
                    terms{i} = [slots{i} ':' vars{i}];
                end
                s = strjoin(terms, ' ');
            else
                s = '';
            end            
        end
        
    end
    
                 
end


