function b = svm_dual_offset(S, K, alpha)
%SVM_DUAL_OFFSET Solve the offset value for dual SVM training
%
%   b = SVM_DUAL_OFFSET(S, K, alpha);
%
%       solves the offset value b for dual SVM training. 
%
%       Note that b is fixed to 0 for SVM ranking.
%
%       Inputs
%       - S:        The SVM problem struct
%       - K:        The kernel matrix
%       - alpha:    The solved alpha vector.
%


% Created by Dahua Lin, on Jan 18, 2012
%

%% main

switch S.type
    case 'class'
        b = offset_for_class(S, K, alpha);
        
    case 'rank'
        b = 0;
        
    otherwise
        error('svm_dual_offset:invalidarg', ...
            'Unsupported SVM type %s for svm_dual_offset.', S.type);
end


%% Core functions

function b = offset_for_class(S, K, alpha)

% select samples on the margin

y = S.y.';
c = S.C.';
ep = 1e-8 * c;
si = find(alpha > ep & alpha < c - ep);

if ~isempty(si)
    sK = K(si, :);
    b = median(y(si) - sK * (y .* alpha));
else
    si = find(alpha > 1e-12);
    if isempty(si)
        b = 0;
        warning('svm_dual_offset:nosuppvec', 'No support vector found.');
        return;
    end
    sy = y(si);
    sV = K(si, :) * (y .* alpha);    
    b0 = median(sy - sV);
    u = sy .* sV;
    
    if isscalar(c)
        b = fminsearch(@(x) fun_c(x, sy, u, c), b0);
    else
        b = fminsearch(@(x) fun_c(x, sy, u, c(si)), b0);
    end
end

function v = fun_c(b, y, u, c)

r = max(1 - u - b * y, 0);
if isscalar(c)
    v = c * sum(r);
else
    v = sum(c .* r);
end
    







    

