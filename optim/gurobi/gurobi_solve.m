function [x, fval, flag, info] = gurobi_solve(P, params)
% Solves a linear or quadratic programming problem using Gurobi solver
%
%   x = gurobi_solve(P);
%   [x, fval] = gurobi_solve(P);
%   [x, fval, flag] = gurobi_solve(P);
%   [x, fval, flag, info] = gurobi_solve(P);
%
%       solves the linear programming problem given by P, where P
%       is a struct representing the problem that can be obtained by
%       calling lp_problem (for LP) or qp_problem (for QP).
%
%       Output arguments:
%       - x:        the solution vector
%       - fval:     the objective value
%       - flag:     the flag indicating the cause of termination
%                   - 0:   solved to optimality
%                   - 1:   proven to be infeasible or unbounded
%                   - 2:   early termination without getting to optima    
%                          
%       - info:     a struct that contains the procedural information, 
%                   which has the following fields:
%                   - status:       the Gurobi status code at termination
%                   - runtime:      the elapsed time in solving the problem
%                   - iters:        the number of simplex iterations
%                   - bariters:     the number of barrier iterations
%                   - objval:       the objective value
%                   - objbound:     the best bound on objective
%
%   ... = gurobi_solve(P, params);
%       
%       additionally specifies the parameters for controlling the solver.
%
%       Here, params should be a struct yielded by the function 
%       gurobi_params.
%

%   History
%   -------
%       - Created by Dahua Lin, on Jan 21, 2011
%


%% verify input

pty = 0;
if isstruct(P) && isfield(P, 'type')
    if strcmp(P.type, 'lp')
        pty = 1;
    elseif strcmp(P.type, 'qp')
        pty = 2;
    end
end
if pty == 0
    error('gurobi_solve:invalidarg', ...
        'The 1st input should be a struct representing an LP or QP problem.');
end

if nargin < 2
    params = [];
else
    if ~(isstruct(params) && numel(params) == 1)
        error('gurobi_solve:invalidarg', ...
            'The params should be a struct scalar.');
    end
end



%% main

% convert to a struct used by the mex core

S.d = double(P.d);
S.f = dvec(P.f);

S.A = dsmat(P.A);
S.bl = dvec(P.bl);
S.bu = dvec(P.bu);

S.Ae = dsmat(P.Aeq);
S.be = dvec(P.beq);

S.lb = dvec(P.l);
S.ub = dvec(P.u);

if pty == 1
    S.Q = [];
else
    S.Q = qmat(P.H);
end

% solve

if nargout <= 1
    x = gurobi_solve_mex(S, params);
elseif nargout == 2
    [x, fval] = gurobi_solve_mex(S, params);
elseif nargout == 3
    [x, fval, flag] = gurobi_solve_mex(S, params);
else
    [x, fval, flag, info] = gurobi_solve_mex(S, params);
end



%% Auxiliary functions

function v = dvec(v)

if ~isempty(v)
    if ~isa(v, 'double') 
        v = double(v);
        if issparse(v)
            v = full(v);
        end
    end
end


function AS = dsmat(A)

if ~isempty(A)
    AS.m = size(A, 1);
    [i, j, v] = find(A);
    AS.i = int32(i) - 1;
    AS.j = int32(j) - 1;
    if ~isa(v, 'double'); v = double(v); end
    AS.v = v;
else
    AS = [];
end
    

function Q = qmat(H)

[I, J, V] = find(H);
Q.I = int32(I) - 1;
Q.J = int32(J) - 1;
if ~isa(V, 'double'); V = double(V); end
Q.V = 0.5 * V;


