function R = smi_varinfer(prg, obs, S0, opts)
% Drive a MCMC sampling procedure based on an SMI program
%
%   R = smi_varinfer(prg, obs, S0, opts);
%
%       Drives a variational inference procedure and collects samples.
%
%       Input arguments:
%       - prg:      the SMI program object, which should be an object
%                   of a class derived from smi_prg.
%
%       - obs:      The observations to be passed to the program.
%
%       - S0:       The initial state to be passed to the program.
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
%   Remarks
%   -------
%       - If opts.nrepeats > 0, multiple optimization procedures will be
%         launched independently. If they use random initialization, or
%         other randomized steps, they would yield different results.
%         Finally, the best one (the one with highest final objective)
%         will be picked and returned.
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 4, 2011
%

%% verify input arguments

if ~isa(prg, 'smi_prg')
    error('smi_varinfer:invalidarg', 'prg should be a SMI program.');
end

opts = varinfer_options(opts);

%% main

nr = opts.nrepeats;

if nr == 1
    R = do_varinfer(prg, obs, S0, opts, 0);
else
    Rs = cell(1, nr);
    parfor i = 1 : nr
        Rs{i} = do_varinfer(prg, obs, S0, opts, i);
    end
    
    % choose the best one
    best_i = 0;
    best_v = -inf;
    
    for i = 1 : nr
        v = Rs{i}.objv(end);
        if v > best_v
            best_i = i;
            best_v = v;
        end            
    end
    R = Rs{best_i};    
end


%% core MCMC function

function R = do_varinfer(prg, obs, S0, opts, ithread)

displevel = opts.display;

% initialization

if displevel >= 2
    fprintf('[[thread %d]] initializing ...\n', ithread);
end

[Sd, Sc] = prg.initialize(obs, S0, 'varinfer');
    
% main iterations

if displevel >= 2
    fprintf('[[thread %d]] iteratively optimizing ...\n', ithread);
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
            fprintf('[[thread %d]]     iter %d\n', ithread, it);
        end
        [Sd, Sc] = prg.update(Sd, Sc);
        
    else
        it_first = it + 1;
        it_last = min(it + ipe, maxiters);
        
        for it = it_first : it_last
            if displevel >= 4
                fprintf('[[thread %d]]     iter %d\n', ithread, it);
            end
            [Sd, Sc] = prg.update(Sd, Sc);
        end
    end
       
    % evaluate objective & determine convergence
    
    cobjv = prg.evaluate_objective(Sd, Sc);
    
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
        fprintf('[[thread %d]]   eval (@ iter %d): objv = %.4g (ch = %g)\n', ...
            ithread, it, cobjv, ch);
    end
            
end

% make output struct

R.sol = prg.make_output(Sd, Sc);
R.niters = it;
R.converged = converged;
R.objv = objv(1:ne);
R.eiters = eiters(1:ne);

