function W = test_mrf2d(m, n, kernel, roi)

if nargin < 4
    roi = empty;
end

use_roi = ~isempty(roi);

N0 = m * n;
h = size(kernel, 1) - 1;
w = size(kernel, 2) - 1;

if use_roi
    N = nnz(roi);
    imap = zeros(m, n);
    imap(roi) = 1 : N;
else
    N = N0;
end

W = zeros(N, N);
for i = 1 : N0
    for j = 1 : N0
        
        if use_roi && ~(roi(i) && roi(j))
            continue;
        end
        
                        
        [iy, ix] = ind2sub([m, n], i);
        [jy, jx] = ind2sub([m, n], j);
        
        dx = abs(jx - ix);
        dy = abs(jy - iy);
                        
        if dx <= w && dy <= h
            if use_roi
                W(imap(i), imap(j)) = kernel(dy+1, dx+1);
            else
                W(i, j) = kernel(dy+1, dx+1);
            end
        end
        
    end
end