classdef smi_state
    %SMI_STATE the base class of the run-time state of an iterative program
    %
    
    % Created by Dahua Lin, on Dec 26, 2011
    %
    
    methods(Abstract)
        
        obj = initialize(obj, varargin);
        % Initialize the state
        
        obj = update(obj);
        % Update current state
        
        R = output(obj);
        % Make the output from the current state        
    
        b = is_ready(obj);
        % Whether the object is properly initialized 
    end
    
end
