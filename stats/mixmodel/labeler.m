classdef labeler
    % The abstract base class for labeling (generating labeling problems)
    %
    
    % Created by Dahua Lin, on April 23, 2010
    %
    
    methods
        
        p = estimate_param(obj, Q)
        % Estimates the labeling parameter from soft assignment                
        
        Pb = accept_param(obj, p);
        % Accept a labeling parameter p and form a labeling problem
        
        L = eval_logpri(obj, p);
        % Estimates the log-prior of the parameter
        
    end
end