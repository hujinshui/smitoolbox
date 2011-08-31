classdef smi_program
    % The class to capture an statistical modeling and inference program
    %
    %   Each program is a collection of compiled code blocks, which 
    %   can refer to each other. A program can be run by an smi driver.
    %
    %   Each block is a struct with the following fields:
    %      
    %   - name:     the name of the block
    %   - steps:    the sequence of steps to execute, which is a struct
    %               array, with following fields:
    %               - op:       the char indicating operation types
    %               - id:       an integer reference ids
    %               - nr:       an repeating-time number
    %               - ivs:      input variable ids
    %               - ovs:      output variable ids
    %
    %
        
    % Created by Dahua Lin, on Aug 29, 2011
    %
    
    %% properties
    
    properties(GetAccess='public', SetAccess='private')
        vars;           % the table of variables
        funcs;          % the table of functions
        blocks;         % the table of blocks
        
        var_map;        % variable name -> id
        func_map;       % function name -> id
        block_map;      % block name -> id
    end
    
    properties(Dependent)
        num_vars;
        num_funcs;
        num_blocks;
    end
    
    methods
        function n = get.num_vars(prg)
            n = numel(prg.vars);
        end
        
        function n = get.num_funcs(prg)
            n = numel(prg.funcs);
        end
        
        function n = get.num_blocks(prg)
            n = numel(prg.blocks);
        end
    end    
    
    
    %% Compilation
    
    methods(Static)
        
        function prg = compile(frmwork)
            % Compiles an SMI framework into a program
            %
            %   prg = smi_program.compile(frmwork);
            %
            %       This function constructs an SMI program by compiling
            %       an SMI framework.
            %
            
            % compile vars
                        
            nv = frmwork.num_vars;
            V = frmwork.vars;
            
            vs = smi_allocate_struct({'name', 'type', 'size'}, nv, 1);
            for i = 1 : nv
                v = V{i};
                vs(i).name = v.name;
                vs(i).type = v.type;
                vs(i).size = v.size;
            end
            vmap = smi_make_map(vs);
            
            % compile functions
            
            nf = frmwork.num_funcs;
            F = frmwork.funcs;
            
            fs = smi_allocate_struct( ...
                {'name', 'obj', 'n_in', 'n_out'}, ...
                nf, 1); 
            
            for i = 1 : nf
                f = F{i};
                fs(i).name = f.name;
                fs(i).obj = f.obj;
                
                Si = f.obj.input_slots;
                So = f.obj.output_slots;
                
                fs(i).n_in = numel(Si);
                fs(i).n_out = numel(So);               
            end
            fmap = smi_make_map(fs);
            
            % compile blocks
            
            B = frmwork.blocks;            
            nb = frmwork.num_blocks;            
            blk_map = smi_make_map(B(1:nb));
                        
            blks = smi_allocate_struct({'name', 'steps'}, nb, 1);
            
            for i = 1 : nb
                blks(i) = smi_program.compile_block( ...
                    B{i}, vs, fs, vmap, fmap, blk_map);
            end
            
            prg = smi_program;
            
            prg.vars = vs;
            prg.funcs = fs;
            prg.blocks = blks;
            
            prg.var_map = vmap;
            prg.func_map = fmap;
            prg.block_map = blk_map;
        end
        
    end
    
       
    methods(Static, Access='private')
        
        function cblk = compile_block(blk, vars, funcs, varmap, funmap, blkmap)
            % Compiles one block
            
            cblk.name = blk.name;
            
            n = blk.num_steps;            
            cblk.steps = smi_allocate_struct( ...
                {'op', 'id', 'nr', 'ivs', 'ovs'}, n, 1);            
            
            for i = 1 : n     
                s = blk.steps{i};
                
                if strcmp(s.type, 'f')
                    
                    cblk.steps(i).op = 'f';
                    
                    f_id = smi_program.get_entity_id( ...
                        funmap, s.func_name, 'function name');
                    cblk.steps(i).id = f_id;
                    
                    cblk.steps(i).nr = s.nrepeats;
                    
                    fun = funcs(f_id);
                    if ~isempty(s.input_slots)
                        ivs = smi_program.get_var_ids( ...
                            fun, 1, vars, varmap, s.input_slots, s.input_vars);
                    else
                        ivs = [];
                    end
                    if ~isempty(s.output_slots)
                        ovs = smi_program.get_var_ids( ...
                            fun, 0, vars, varmap, s.output_slots, s.output_vars);
                    else
                        ovs = [];
                    end
                    
                    if ~fun.obj.test_slots(ivs > 0, ovs > 0)
                        error('smi_program:compile_error', ...
                            ['In step %d of block %s, ', ...
                            'the function is not properly connected'], ...
                            i, blk.name);
                    end
                    
                    cblk.steps(i).ivs = ivs;
                    cblk.steps(i).ovs = ovs;              
                
                elseif strcmp(s.type, 'b')
                    
                    cblk.steps(i).op = 'b';
                    
                    b_id = smi_program.get_entity_id( ...
                        blkmap, s.block_name, 'block name');
                    cblk.steps(i).id = b_id;
                    
                    cblk.steps(i).nr = s.nrepeats;
                    cblk.steps(i).ivs = [];
                    cblk.steps(i).ovs = [];
                                        
                else
                    error('smi_program:compile_error', ...
                        'Invalid step type %s', s.type);
                end
            end                                    
        end
        
    end
    
    
    
    %% Displaying functions
    
    methods
        
        function disp(prg)
            % Display the program
            %
            %   disp(prg);
            %
            
            fprintf('SMI Program\n');
            fprintf('------------------\n');
            
            fprintf('  # Variables = %d\n', prg.num_vars);
            fprintf('  # Functions = %d\n', prg.num_funcs);
            fprintf('  # Blocks = %d\n', prg.num_blocks);
            fprintf('\n');
            
        end
        
        
        function dump(prg, fid)
            % Display detailed information about the program
            %
            %   dump(prg);
            %   dump(prg, fid);
            %
            
            if nargin < 2
                fid = 1;
            end
            
            fprintf(fid, 'SMI Program (Details)\n');
            fprintf(fid, '===============================\n');
            
            fprintf(fid, 'Variables\n');
            fprintf(fid, '----------------\n');
            
            nv = prg.num_vars;
            V = prg.vars;
            for i = 1 : nv
                v = V(i);
                fprintf(fid, '    [%d] %s: %s %s\n', ...
                    i, v.name, v.type, smi_vsize2str(v.size));
            end
            fprintf(fid, '\n');
            
            fprintf(fid, 'Functions\n');
            fprintf(fid, '----------------\n');
            
            nf = prg.num_funcs;
            F = prg.funcs;
            for i = 1 : nf
                f = F(i);
                fprintf(fid, '    [%d] %s: class %s (%d in, %d out)\n', ...
                    i, f.name, class(f.obj), f.n_in, f.n_out);
                
                inlets = f.obj.input_slots;
                outlets = f.obj.output_slots;
                
                for j = 1 : numel(inlets)
                    sl = inlets(j);
                    fprintf(fid, '%s  <in> %-15s: %s %s\n', ...
                        blanks(8), sl.name, sl.type, smi_vsize2str(sl.size));
                end

                for j = 1 : numel(outlets)
                    sl = outlets(j);
                    fprintf(fid, '%s <out> %-15s: %s %s\n', ...
                        blanks(8), sl.name, sl.type, smi_vsize2str(sl.size));
                end
                
                fprintf(fid, '\n');
            end
            
            fprintf(fid, 'Blocks\n');
            fprintf(fid, '----------------\n');
            
            nb = prg.num_blocks;
            B = prg.blocks;
            
            for i = 1 : nb
                blk = B(i);
                ns = numel(blk.steps);
                fprintf(fid, 'Block %s (len = %d)\n', blk.name, ns);
                
                for j = 1 : ns
                    
                    st = blk.steps(j);                   
                    
                    if st.op == 'f'
                        
                        fname = F(st.id).name;
                        fobj = F(st.id).obj;
                        
                        islots = fobj.input_slots;
                        oslots = fobj.output_slots;
                        
                        n_in = F(st.id).n_in;
                        n_out = F(st.id).n_out;
                        
                        fprintf(fid, '    [%d] call func %s [#=%d]: (%s) => (%s)\n', ...
                            j, fname, st.nr, ...
                            smi_program.make_connstr(n_in, islots, V, st.ivs), ...
                            smi_program.make_connstr(n_out, oslots, V, st.ovs));                         
                        
                    elseif st.op == 'b'
                        
                        bname = B(st.id).name;
                        
                        fprintf(fid, '    [%d] call block %s [#=%d]\n', ...
                            j, bname, st.nr);                        
                        
                    end                    
                end                
                fprintf('\n');
            end                        
        end    
        
    end 
    
    
    
    methods(Static, Access='private')
        
        function s = make_connstr(n, slots, vars, vids)
            % make a string to represent connections
            
            if any(vids > 0)                
                terms = cell(1, n);
                for i = 1 : n
                    vid = vids(i);
                    if vid > 0                     
                        terms{i} = [slots(i).name ':' vars(vid).name];
                    end                
                end            
                s = strjoin(terms(vids > 0), ', ');
            else
                s = '';
            end
            
        end
    
    end
    
    
    
    %% Auxiliary functions
       
    methods(Static, Access='private')                
        
        function id = get_entity_id(M, name, title)
            % Gets the id of an entity from a map
            
            if isfield(M, name)
                id = M.(name);
            else
                error('There is no %s %s in the program.', title, name);
            end
        end
        
        function vids = get_var_ids(fun, is_input, vars, vmap, cslots, cvars)
            
            % Extract the ids to form an array of vars

            if is_input
                fslots = fun.obj.input_slots;
            else
                fslots = fun.obj.output_slots;
            end
            vids = zeros(1, numel(fslots));
            slotnames = {fslots.name};            
            
            for i = 1 : numel(cslots)
                sn = cslots{i};
                vn = cvars{i};
                
                [tf, islot] = ismember(sn, slotnames);
                if ~tf
                    error('There is no slot in function %s called %s', ...
                        funname, s);
                end
                slot = fslots(islot);
                
                tf = isfield(vmap, vn);
                if ~tf
                    error('There is no variable name %s in the program', v);
                end
                ivar = vmap.(vn);
                va = vars(ivar);
               
                if vids(islot) ~= 0
                    error('The slot %s of function %s is connected more than once.', ...
                        s, fun.name);
                end
                
                if is_input
                    if ~smi_is_vcompatible(va, slot)
                        error('Incompatible conversion: %s => %s\n', ...
                            ['variable ', va.name], ...
                            ['slot ', slot.name, ' of function ', fun.name]);
                    end
                else
                    if ~smi_is_vcompatible(slot, va)
                        error('Incompatible conversion: %s <= %s\n', ...
                            ['variable ', va.name], ...
                            ['slot ', slot.name, ' of function ', fun.name]);
                    end
                end
                
                vids(islot) = ivar;
            end
        end
        
        
    end
    
    
end
