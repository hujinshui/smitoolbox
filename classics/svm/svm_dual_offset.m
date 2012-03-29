function b = svm_dual_offset(S, K, alpha, si)
%SVM_DUAL_OFFSET Solve the offset value for dual SVM training
%
%   b = SVM_DUAL_OFFSET(S, K, alpha, si);
%
%       solves the offset value b for dual SVM training. 
%
%       Note that b is fixed to 0 for SVM ranking.
%
%       Inputs
%       - S:        The SVM problem struct
%       - K:        The kernel matrix over the entire training set
%       - alpha:    The solved alpha vector.
%       - si:       The indices of support vectors.
%


% Created by Dahua Lin, on Jan 18, 2012
%

%% main

switch S.type
    case 'class'
        b = offset_for_class(S, K, alpha, si);
        
    case 'regress'
        b = offset_for_regress(S, K, alpha, si);
        
    case 'rank'
        b = 0;
        
    otherwise
        error('svm_dual_offset:invalidarg', ...
            'Unsupported SVM type %s for svm_dual_offset.', S.type);
end


%% Core functions

function b = offset_for_class(S, K, alpha, si)

y = S.y.';
c = S.C.';
ep = min(1e-9 * c, 1e-6);
ui = find(alpha > ep & alpha < c - ep);
a = y(si) .* alpha(si);

if ~isempty(ui)
    b = median(y(ui) - K(ui, si) * a);
else
    b = 0;
    warning('svm_dual_offset:nosuppvec', ...
        'No training samples found residing on the margin.');
end
    

function b = offset_for_regress(S, K, alpha, si)

n = S.n;
a1 = alpha(1:n);
a2 = alpha(n+1:end);

if isscalar(S.C)
    c = S.C;
else
    c = S.C.';
end

ep = min(1e-9 * c, 1e-6);
s1 = find(a1 < c - ep | a2 > ep);
s2 = find(a1 > ep | a2 < c - ep);

if ~isempty(s1) && ~isempty(s2)
    tol = S.tol;
    y = S.y.';    
    ad = a1 - a2;
    
    Ks = K(:, si);    
    ad = ad(si);  
    
    m1 = max(y(s1) - Ks(s1, :) * ad - tol);
    m2 = min(y(s2) - Ks(s2, :) * ad - tol);    
    b0 = (m1 + m2) / 2;   
           
    u = Ks * ad;    
    b = fminsearch(@(x) fun_r(x, u, y, c), b0);    
else
    b = 0;
    warning('svm_dual_offset:nosuppvec', ...
        'No proper training samples to bound of value of b.');
end


function v = fun_r(b, u, y, c)
% objective fun for searching b

loss = max(abs(u + b - y), 0);
if isscalar(c)
    v = c * sum(loss);
else
    v = sum(loss .* c);
end

