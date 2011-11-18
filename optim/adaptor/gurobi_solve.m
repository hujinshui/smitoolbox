function [x, fval, flag, info] = gurobi_solve(P, options)
% Adaptor of GUROBI optimizer to solve the LP/QP problem 
%
%   x = gurobi_solve(P);
%   x = gurobi_solve(P, options);
%   
%       Uses GUROBI optimizer to solve the LP or QP problem given by the 
%       struct P (as constructed in lp_problem or qp_problem).
%
%   [x, fval] = gurobi_solve( ... );
%   [x, fval, flag] = gurobi_solve( ... );
%   [x, fval, flag, info] = gurobi_solve( ... );
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
    error('gurobi_solve:invalidarg', 'P should be a qp_problem struct.');
end

is_qp = strcmp(P.type, 'qp');

if nargin < 2
    options = [];
else
    if ~isstruct(options)
        error('gurobi_solve:invalidarg', 'options should be a GUROBI params struct.');
    end
end

if exist('gurobi_lpqp', 'file') ~= 2
    error('gurobi_solve:invalidarg', ...
        'This function relies on gurobi_lpqp (a gurobi matlab wrapper), which is not found.');
end


%% main

% convert problem

if is_qp
    GP.H = P.H;
end

GP.f = P.f;

GP.A = P.A;
GP.bl = P.bl;
GP.bu = P.bu;

GP.Ae = P.Aeq;
GP.be = P.beq;

GP.lb = P.l;
GP.ub = P.u;

% solve

nout = nargout;

if nout <= 1
    x = gurobi_lpqp(GP, options);
elseif nout == 2
    [x, fval] = gurobi_lpqp(GP, options);
elseif nout == 3
    [x, fval, flag] = gurobi_lpqp(GP, options);
else
    [x, fval, flag, info] = gurobi_lpqp(GP, options);
end




