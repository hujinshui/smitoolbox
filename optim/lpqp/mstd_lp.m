function [x, fval, flag, info] = mstd_lp(P, options)
% Adaptor of linprog to solve the LP problem given by lp_problem struct
%
%   x = mstd_lp(P);
%   x = mstd_lp(P, options);
%   
%       Uses linprog (in optimization toolbox) to solve the LP problem
%       given by the struct P (as constructed in lp_problem).
%
%   [x, fval] = mstd_lp( ... );
%   [x, fval, flag] = mstd_lp( ... );
%   [x, fval, flag, info] = mstd_lp( ... );
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

if ~(isstruct(P) && isfield(P, 'type') && strcmp(P.type, 'lp'))
    error('mstd_lp:invalidarg', 'P should be a lp_problem struct.');
end

if nargin < 2
    options = optimset('Display', 'off');
else
    if ~isstruct(options)
        error('mstd_lp:invalidarg', 'options should be a struct.');
    end
end

%% main

% convert problem

f = P.f;

if isempty(P.A)
    A = [];
    b = [];
else
    if isempty(bl)
        A = P.A;
        b = P.bu;
    elseif isempty(bu)
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

if nout <= 1
    x = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
elseif nout == 2
    [x, fval] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
elseif nout == 3
    [x, fval, flag] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
else
    [x, fval, flag, info] = linprog(f, A, b, Aeq, beq, lb, ub, [], options);
end


    





