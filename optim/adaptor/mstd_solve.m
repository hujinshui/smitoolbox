function [x, fval, flag, info] = mstd_solve(P, options)
% Adaptor of MATLAB optimization toolbox to solve the LP/QP problem
%
%   x = mstd_solve(P);
%   x = mstd_solve(P, options);
%   
%       Uses linprog/quadprog (in optimization toolbox) to solve the 
%       LP/QP problem given by the struct P (as constructed in lp_problem 
%       or qp_problem).
%
%   [x, fval] = mstd_solve( ... );
%   [x, fval, flag] = mstd_solve( ... );
%   [x, fval, flag, info] = mstd_solve( ... );
%
%       Returns additional outputs.
%
%       - fval: the value of objective function at x
%       - flag: the exit flag that describes the exit condition.
%       - info: the struct that contains information of the optimization
%               procedure.
%
%       Please refer to the document of linprog for more details about
%       the options and outputs.
%

% Created by Dahua Lin, on Apr 7, 2011
%

%% verify input arguments

if ~(isstruct(P) && isfield(P, 'type') && (strcmp(P.type, 'lp') || strcmp(P.type, 'qp')))
    error('mstd_solve:invalidarg', 'P should be a lp_problem/qp_problem struct.');
end

is_qp = strcmp(P.type, 'qp');

if nargin < 2
    if ~is_qp
        options = optimset('linprog');
    else
        options = optimset('quadprog');
    end
    options = optimset(options, 'Display', 'off');
else
    if ~isstruct(options)
        error('mstd_solve:invalidarg', 'options should be a struct.');
    end
end

%% main

% convert problem

f = P.f;

if is_qp
    H = P.H;
end

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

nout = nargout;

if ~is_qp
    if nout <= 1
        x = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 2
        [x, fval] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 3
        [x, fval, flag] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
    else
        [x, fval, flag, info] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
    end
else
    if nout <= 1
        x = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 2
        [x, fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    elseif nout == 3
        [x, fval, flag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    else
        [x, fval, flag, info] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
    end
end


    





