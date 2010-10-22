function W = test_mrf2d(siz, kernel)


m = siz(1);
n = siz(2);

h = size(kernel, 1) - 1;
w = size(kernel, 2) - 1;

N = m * n;


W = zeros(N, N);
for i = 1 : N
    for j = 1 : N
                                
        [iy, ix] = ind2sub([m, n], i);
        [jy, jx] = ind2sub([m, n], j);
        
        dx = abs(jx - ix);
        dy = abs(jy - iy);
                        
        if dx <= w && dy <= h
            W(i, j) = kernel(dy+1, dx+1);
        end
        
    end
end