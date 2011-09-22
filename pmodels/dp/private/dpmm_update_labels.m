function [A, L, C, Z, K, ch] = dpmm_update_labels(mdl, alpha, K, A, L, L0, priC, C, Z, inds)
% The implementation of DPMM label updating function
%
%   arguments:
%       - mdl:      the model (of class nonparam_model)
%       - alpha:    the concentration parameter
%       - L:        the log likelihood
%       - L0:       the log likelihood with respect to base
%       - priC:     the prior counts (can be [])
%       - A:        the cell array of atoms
%       - C:        the current counts
%       - Z:        labels
%       - inds:     the sequence of indices to be updated
%

%% main

len = numel(inds);
if len == 0
    return;
end

Kpri = numel(priC);

% pre-computation

pri_ = C.';
if Kpri > 0
    pri_ = pri_(1:Kpri) + priC.';
end
lpri = log(pri_);

lalpha = log(alpha);

rnums = rand(1, len);
ch = false;

for j = 1 : len    
    i = inds(j);
    
    % remove old label
    z0 = Z(i);
    if z0 > 0
        Z(i) = 0;
        C(z0) = C(z0) - 1;
    end
    
    % draw new label
    if K == 0
        z = 1;
    else
        ev = lpri(1:K) + L(1:K, i);
        e0 = lalpha + L0(i);
        
        mv = max(max(ev), e0);
        
        w = [exp(ev - mv); exp(e0 - mv)];
        sw = sum(w);
                
        tv = rnums(j) * sw;
        cw = 0;
        z = 0;
        
        while cw < tv
            z = z + 1;
            cw = cw + w(z);
        end        
    end
    
    % add new atoms if necessary
    
    if z > K  % z = K + 1
        a = mdl.create_atom(i);
        
        capa = numel(A);
        if capa == K
            A{1, 2 * K} = [];
            C(1, 2 * K) = 0;
            L(2 * K, end) = 0;
            
            lpri(2 * K, 1) = 0;
        end
        
        K = K + 1;
        A{K} = a;
        L(K, :) = mdl.evaluate_loglik(a);        
    end
    
    % set label
    
    Z(i) = z;
    C(z) = C(z) + 1;    
    
    if z0 ~= z
        if z0 >= 1
            if z0 > Kpri
                lpri(z0) = log(C(z0));
            else
                lpri(z0) = log(C(z0) + priC(z0));
            end
        end
        if z > Kpri
            lpri(z) = log(C(z));
        else
            lpri(z) = log(C(z) + priC(z));
        end
        
        ch = true;
    end
    
end

