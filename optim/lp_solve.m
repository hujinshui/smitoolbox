function [x, info] = lp_solve(prob, solver, optimParam, printLevel)
% Solve a Linear Programming problem using specified solver and parameters
%
%   x = lp_solve(prob, solver);
%   x = lp_solve(prob, solver, optimParam);
%   x = lp_solve(prob, solver, optimParam, printLevel);
%
%       solves a linear programming problem prob using the specified
%       solver. 
%
%       optimParam is a struct whose fields specifies the optimization
%       options. printLevel is an integer value indicating the level of
%       details in printing.
%
%       Here is a list of supported solvers:
%
%       Solvers in MATLAB Optimization toolbox:
%
%       'matlab.activeset': The active-set algorithm in linprog
%       'matlab.simplex':   The simplex algorithm 
%       'matlab.lipsol':    The linear interior point algorithm (LIPSOL)
%
%       Solver in MOSEK:
%
%       'mosek.default':    The default solver in MOSEK
%
%       Solvers in TOMLAB that can solve LP problems:
%
%       'tomlab.lpsimplex'
%       'tomlab.minos'     
%       'tomlab.lpopt'
%       'tomlab.cplex'
%       'tomlab.gurobi'
%

%   History
%   -------
%       - Created by Dahua Lin, on Sep 25, 2010
%

%% verify and parse input

if ~(isstruct(prob) && isfield(prob, 'type') && strcmp(prob.type, 'lp'))
    error('lp_solve:invalidarg', ...
        'prob should be a LP problem struct.');
end

if ~ischar(solver)
    error('lp_solve:invalidarg', 'solver should be a string.');
end

[pkgname, solvername] = parse_solver_name(lower(solver));

if nargin >= 3 && ~isempty(optimParam)
    if ~isstruct(optimParam)
        error('lp_solve:invalidarg', 'optimParam should be a struct.');
    end
else
    optimParam = [];
end

if nargin >= 4 && ~isempty(printLevel)
    if ~(isnumeric(printLevel) && isscalar(printLevel) && printLevel >= 0)
        error('lp_solve:invalidarg', 'printLevel should be a non-negative integer.');
    end
else
    printLevel = 0;
end

%% main

switch pkgname
    case 'matlab'
        if nargout < 2
            x = lp_solve_matlab(prob, solvername, optimParam, printLevel);
        else
            [x, info] = lp_solve_matlab(prob, solvername, optimParam, printLevel);
        end
    case 'mosek'
        if nargout < 2
            x = lp_solve_mosek(prob, solvername, optimParam, printLevel);
        else
            [x, info] = lp_solve_mosek(prob, solvername, optimParam, printLevel);
        end
    case 'tomlab'
        if nargout < 2
            x = lp_solve_tomlab(prob, solvername, optimParam, printLevel);
        else
            [x, info] = lp_solve_tomlab(prob, solvername, optimParam, printLevel);
        end
    otherwise
        error('lp_solve:invalidarg', 'Invalid solver name %s', solver);
end
    

%% core functions

function [x, info] = lp_solve_matlab(prob, solvername, optimParam, printLevel)

[A, b, Aeq, beq] = to_matlab_constrs(prob);
dispstr = to_matlab_dispstr(printLevel);

switch solvername
    case 'activeset'
        largeScale = 'off';
        simplex = 'off';
    case 'simplex'
        largeScale = 'off';
        simplex = 'on';
    case 'lipsol'
        largeScale = 'on';
        simplex = 'off';        
    otherwise
        error('lp_solve:invalidarg', ...
            'Invalid solver name %s', ['matlab.' solvername]);
end

if isempty(optimParam)
    optimParam = optimset('linprog');
end
optimParam = optimset(optimParam, ...
    'Display', dispstr, ...
    'LargeScale', largeScale, ...
    'Simplex', simplex);

if nargout < 2
    x = linprog(prob.c, A, b, Aeq, beq, prob.l, prob.u, prob.x0, optimParam);
else
    [x, fval, exitflag, info] = ...
        linprog(prob.c, A, b, Aeq, beq, prob.l, prob.u, prob.x0, optimParam);
    
    info.objval = fval;
    info.exitflag = exitflag;
    info.converged = exitflag > 0;
end


function [x, info] = lp_solve_mosek(prob, solvername, optimParam, printLevel)

if ~strcmp(solvername, 'default')
    error('lp_solve:invalidarg', ...
        'Invalid solver name %s', ['mosek.' solvername]);
end

[A, bl, bu] = to_tomlab_constrs(prob);

mp.c = prob.c;
if issparse(A)
    mp.a = A;
else
    mp.a = sparse(A);
end
mp.blc = bl;
mp.buc = bu;
mp.blx = prob.l;
mp.bux = prob.u;

optimParam.Display = to_mosek_dispstr(printLevel);

[cmd, verb, param] = msksetup(1, optimParam); %#ok<ASGLU>
[rcode, info] = mosekopt(cmd, mp, param); %#ok<ASGLU>

if isfield(info, 'sol')
  x = info.sol.itr.xx;
else
  x = [];
end


function [x, info] = lp_solve_tomlab(prob, solvername, optimParam, printLevel)

[A, bl, bu] = to_tomlab_constrs(prob);

% construct problem

tp = lpAssign(prob.c, A, bl, bu, prob.l, prob.u, prob.x0);
tp.PriLevOpt = printLevel;
if ~isempty(optimParam)
    tp.optParam = optimParam;
end

% solve

switch solvername
    case 'lpsimplex'
        info = lpSimplex(tp);
        x = info.x_k;
                
    case 'minos'
        tp = ProbCheck(tp,'minos',8);
        info = minoslpTL(tp);
        x = info.x_k;
        
    case 'lpopt'
        tp = ProbCheck(tp,'lpopt',8);
        info = lpoptTL(tp);
        x = info.x_k;
        
    case 'cplex'
        tp = ProbCheck(tp,'cplex',8);
        info = cplexTL(tp);
        x = info.x_k;
        
    case 'gurobi'   
        tp = ProbCheck(tp,'gurobi',8);
        info = gurobiTL(tp);
        x = info.x_k;
        
    otherwise
        error('lp_solve:invalidarg', ...
            'Invalid solver name %s', ['tomlab.' solvername]);
end


