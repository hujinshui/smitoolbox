function C = confusmat(m, L1, L2, op)
% Compute the confusion matrix between two sets of labels
%
%   C = confusmat(m, L1, L2);
%       computes the confusion matrix between the labels given in L1
%       and that given in L2.
%
%       L1 and L2 should be vectors of the same size, and m is the
%       number of classes. The value range of labels should be in 
%       [1, m].
%
%       Then C is a matrix of size m x m, where C(k, l) is the number
%       of samples with L1(i) == k and L2(i) == l.
%
%   C = confusmat([m1, m2], L1, L2);
%       computes the confusion matrix between the labels in L1 and L2,
%       where L1 and L2 are respectively based on m1 and m2 classes.
%
%       Actually, confusmat(m, L1, L2) is equivalent to 
%       confusmat([m m], L1, L2).
%
%   C = confusmat(m, L1, L2, 'nrm')
%       computes the normalized confusion matrix, of which the sum of
%       all entries are normalized to one.
%
%       A normalized confusion matrix can be considered as the probability
%       mass function of the joint distribution.
%
%   Remarks
%   -------
%       - The caller should ensure that the labels are actually in the
%         range of [1, m]. Otherwise, the function will raise and error.
%

%   History
%   -------
%       - Created by Dahua Lin, on May 26,2010
%

%% verify input arguments

if ~(isnumeric(m) && (isscalar(m) || numel(m) == 2))
    error('confusmat:invalidarg', 'm should be a either a scalar or a pair.');   
end

if isscalar(m)
    m1 = m;
    m2 = m;
else
    m1 = m(1);
    m2 = m(2);
end

if ~(isvector(L1) && isvector(L2) && isnumeric(L1) && isnumeric(L2) && ...
        isequal(size(L1), size(L2)))
    error('confusmat:invalidarg', ...
        'L1 and L2 should be numeric vectors of the same size.');
end

if any(L1 > m1) || any(L1 < 1)
    error('confusmat:invalidarg', 'Some labels in L1 exceed range.');
end
if any(L2 > m2) || any(L2 < 1)
    error('confusmat:invalidarg', 'Some labels in L2 exceed range.');
end

nrm = false;
if nargin >= 4
    if strcmp(op, 'nrm')
        nrm = true;
    else
        error('confusmat:invalidarg', ...
            'The 4th argument can only be ''nrm''.');
    end
end


%% main

I = L1 + m1 * (L2 - 1);
C = intcount([1, m1 * m2], I);

if nrm
    C = C / sum(C);
end

C = reshape(C, m1, m2);


