function R = varinfer_drive(S, opts)
% Drive a MCMC sampling procedure based on an SMI program
%
%   R = varinfer_drive(S, opts);
%
%       Drives a variational inference procedure and collects samples.
%
%       Input arguments:
%       - S:        the SMI state that has been properly initialized.
%
%       - opts:     the inference control options, which can be obtained 
%                   by calling varinfer_options.
%                   (See the help of varinfer_options for details).
%       
%       Output arguments:
%       - R:        a struct with following fields:
%                   - 'sol':    the obtained solution (this is the
%                               result returned by invoking make_output
%                               method on the states of final iteration)
%
%                   - 'state':  the updated state object
%
%                   - 'niters': the number of elapsed iterations
%
%                   - 'converged': whether the procedure converges
%
%                   - 'objv':   the vector of recorded objective values.
%                               objv(end) is the final objective.
%
%                   - 'eiters': the indices of iterations at which the
%                               objectives are evaluated. Specifically,
%                               objv(i) is evaluated at iteration eiters(i).
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 4, 2011
%       - Modified by Dahua Lin, on Dec 27, 2011
%

%% verify input arguments

if ~isa(S, 'smi_state')
    error('varinfer_drive:invalidarg', 'S should be an smi_state object.');
end

if ~S.is_ready()
    error('varinfer_drive:invalidarg', 'S has not been ready.');
end

opts = varinfer_options(opts);

%% main

displevel = opts.display;
    
% main iterations

if displevel >= 2
    fprintf('Variational inference updating ...\n');
end

converged = false;
it = 0;
maxiters = opts.maxiters;
ipe = opts.ipe;

ne = 0;
max_ne = ceil(maxiters / ipe); 

objv = zeros(1, max_ne);
eiters = zeros(1, max_ne);

while ~converged && it < maxiters
    
    % do update
    
    if ipe == 1        
        it = it + 1;
        if displevel >= 4
            fprintf('     iter %d\n', it);
        end
        S = S.update();
        
    else
        it_first = it + 1;
        it_last = min(it + ipe, maxiters);
        
        for it = it_first : it_last
            if displevel >= 4
                fprintf('     iter %d\n', it);
            end
            S = S.update();
        end
    end
       
    % evaluate objective & determine convergence
    
    cobjv = S.evaluate_objv();
    
    ne = ne + 1;
    eiters(ne) = it;
    objv(ne) = cobjv;
    if ne == 1
        ch = nan;
    else
        ch = objv(ne) - objv(ne-1);
        converged = (abs(ch) <= opts.tol);
    end
    
    if displevel >= 3
        fprintf('   eval (@ iter %d): objv = %.4g (ch = %g)\n', ...
            it, cobjv, ch);
    end
            
end

% make output struct

R.sol = S.output();
R.state = S;
R.niters = it;
R.converged = converged;
R.objv = objv(1:ne);
R.eiters = eiters(1:ne);



