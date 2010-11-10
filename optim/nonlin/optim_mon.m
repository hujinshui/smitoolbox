classdef optim_mon
    % The class to represent a optimization monitor
    %
    
    % Created by Dahua Lin, on Aug 7, 2010
    %
    
    properties(Constant)        
        NoneLevel = 0;
        ProcLevel = 1;
        IterLevel = 2;
    end
        
    properties(GetAccess='public', SetAccess='private')
        level;
        fid = 1;
    end
    
    methods
        
        function mon = optim_mon(level)
            % Constructs an optimization monitor
            %
            %   mon = optim_mon(level);
            %       constructs an optimization monitor at specified level.
            %       Here, level can be either a number of a string name:
            %       'none', 'proc', 'iter', 'step', or 'term'.
            %
            
            if isnumeric(level) && isscalar(level)
                v = level;
            elseif ischar(level)
                switch lower(level)
                    case 'none'
                        v = 0;
                    case 'proc'
                        v = 1;
                    case 'iter'
                        v = 2;
                    case 'step'
                        v = 3;
                    otherwise
                        error('optim_mon:invalidarg', ...
                            'Unknown monitor level %s', level);
                end
            else
                error('optim_mon:invalidarg', ...
                    'The input monitor level is invalid.');
            end
            
            mon.level = v;
        end
    
        
        function on_proc_start(mon)
            if mon.level >= optim_mon.ProcLevel
                if mon.level >= optim_mon.IterLevel
                    print_log(mon, '%10s%15s%17s%15s\n', ...
                        'Iters', 'Move', 'Obj-value', 'Obj-change');
                end
            end
        end        
        
        function on_proc_end(mon, info)
            if mon.level >= optim_mon.ProcLevel
                if info.IsConverged
                    print_log(mon, 'Optimization converges with %d iterations.\n', ...
                        info.NumIters);
                else
                    print_log(mon, ...
                        'Optimization did NOT converge after %d iterations.\n', ...
                        info.NumIters);
                end
            end
        end
        
        function on_iter_start(mon, it) %#ok<MANU,INUSD>            
        end
        
        function on_iter_end(mon, it, itstat)            
            if mon.level >= optim_mon.IterLevel
                print_log(mon, '%10d  %13.4g  %15.6g  %13.4g\n', ...
                    it, itstat.MoveNorm, itstat.FunValue, itstat.FunChange);
            end
        end
                        
    end

    
    methods(Access='private')
        
        function print_log(mon, varargin)
            fprintf(mon.fid, varargin{:});
        end
        
    end
            
end



