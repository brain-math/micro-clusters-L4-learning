function W2 = circmat_align_2D(W, ifForward)
% W: (N, N); N_sqrt should be odd
% ifForward == 1: from circulant to centered
% ifForward == 0: from centered back to circulant

N = size(W, 1); N_sqrt = sqrt(N);
idx_ctr = (N_sqrt + 1) / 2;

W2 = NaN(N, N);
for xi = 1: N
    i = mod(xi - 1, N_sqrt) + 1; j = ceil(xi / N_sqrt);
    Wk_2D = reshape(W(xi, :), [N_sqrt N_sqrt]);
    if ifForward == 1
        W2k_2D = circshift(Wk_2D, [idx_ctr - i, idx_ctr - j]);
    elseif ifForward == 0
        W2k_2D = circshift(Wk_2D, [i - idx_ctr, j - idx_ctr]);
    end
    W2(xi, :) = reshape(W2k_2D, [1 N]);
end
