function [x, fval, flag, info] = cplex_solve(P, options)
% Adaptor of CPLEX optimizer to solve the LP/QP problem 
%
%   x = cplex_solve(P);
%   x = cplex_solve(P, options);
%   
%       Uses quadprog (in optimization toolbox) to solve the LP or QP 
%       problem given by the struct P (as constructed in lp_problem or 
%       qp_problem).
%
%   [x, fval] = cplex_solve( ... );
%   [x, fval, flag] = cplex_solve( ... );
%   [x, fval, flag, info] = cplex_solve( ... );
%
%       Returns additional outputs.
%
%       - fval: the value of objective function at x
%       - flag: the exit flag that describes the exit condition.
%       - info: the struct that contains information of the optimization
%               procedure.
%

% Created by Dahua Lin, on Apr 15, 2011
%

%% verify input arguments

if ~(isstruct(P) && isfield(P, 'type') && (strcmp(P.type, 'lp') || strcmp(P.type, 'qp')))
    error('cplex_solve:invalidarg', 'P should be a qp_problem struct.');
end

is_qp = strcmp(P.type, 'qp');

if nargin < 2
    options = [];
else
    if ~isstruct(options)
        error('cplex_solve:invalidarg', 'options should be a CPLEX option struct.');
    end
end

%% main

% convert problem

H = P.H;
f = P.f;

if isempty(P.A)
    A = [];
    b = [];
else
    if isempty(P.bl)
        A = P.A;
        b = P.bu;
    elseif isempty(P.bu)
        A = -P.A;
        b = -P.bl;
    else
        A = [P.A; -P.A];
        b = [P.bu; -P.bl];
    end
end

Aeq = P.Aeq;
beq = P.beq;

lb = P.l;
ub = P.u;

% solve

if isempty(options);
    options = cplexoptimset('Display', 'off');
end
            
    
nout = nargout;

if ~is_qp
    if nout <= 1
        x = cplexlp(f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 2
        [x, fval] = cplexlp(f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 3
        [x, fval, flag] = cplexlp(f, A, b, Aeq, beq, lb, ub, [], options);
    else
        [x, fval, flag, info] = cplexlp(f, A, b, Aeq, beq, lb, ub, [], options);
    end
else
    if nout <= 1
        x = cplexqp(H, f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 2
        [x, fval] = cplexqp(H, f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 3
        [x, fval, flag] = cplexqp(H, f, A, b, Aeq, beq, lb, ub, [], options);
    else
        [x, fval, flag, info] = cplexqp(H, f, A, b, Aeq, beq, lb, ub, [], options);
    end
end


