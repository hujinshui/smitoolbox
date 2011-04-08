function [x, fval, flag, info] = mstd_qp(P, options)
% Adaptor of quadprog to solve the QP problem given by lp_problem struct
%
%   x = mstd_qp(P);
%   x = mstd_qp(P, options);
%   
%       Uses quadprog (in optimization toolbox) to solve the QP problem
%       given by the struct P (as constructed in qp_problem).
%
%   [x, fval] = mstd_qp( ... );
%   [x, fval, flag] = mstd_qp( ... );
%   [x, fval, flag, info] = mstd_qp( ... );
%
%       Returns additional outputs.
%
%       - fval: the value of objective function at x
%       - flag: the exit flag that describes the exit condition.
%       - info: the struct that contains information of the optimization
%               procedure.
%
%       Please refer to the document of quadprog for more details about
%       the options and outputs.
%

% Created by Dahua Lin, on Apr 7, 2011
%

%% verify input arguments

if ~(isstruct(P) && isfield(P, 'type') && strcmp(P.type, 'qp'))
    error('mstd_qp:invalidarg', 'P should be a qp_problem struct.');
end

if nargin < 2
    options = [];
else
    if ~isstruct(options)
        error('mstd_qp:invalidarg', 'options should be a struct.');
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

if isempty(options);
    if ( ...
            (isempty(A) && isempty(Aeq)) || ...
            (isempty(A) && isempty(lb) && isempty(ub) && size(Aeq,1) <= P.d) )
        
        options = optimset('LargeScale', 'on', 'Display', 'off');
    else
        options = optimset('LargeScale', 'off', 'Display', 'off');
    end
end
            
    
nout = nargout;

if nout <= 1
    x = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
elseif nout == 2
    [x, fval] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
elseif nout == 3
    [x, fval, flag] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
else
    [x, fval, flag, info] = quadprog(H, f, A, b, Aeq, beq, lb, ub, [], options);
end


