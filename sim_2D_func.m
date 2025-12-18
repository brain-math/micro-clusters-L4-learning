function sim_2D_func(N_sqrt, rho_deg, alpha_width_C, alpha_width_K, sigma_A_input_deg,...
    N_step, idx_trace, savename, ifUse73Save, ifDisplay)

%% Coordinates
% N_sqrt = 21;    % must be odd
N = N_sqrt ^ 2;
%
x_sqrt = -(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2;    % unit: grid
[j_grid, i_grid] = meshgrid(x_sqrt);
% used as (x1, x2) or (alpha1, alpha2) -- use ij axis, not xy!
% So that later you can use reshape to convert among 1D, 2D, 4D.
r_grid = sqrt(j_grid .^ 2 + i_grid .^ 2);
idx_ctr = (N_sqrt + 1) / 2;
%
rho_mag = 20;    % (um / deg)    % physiological parameter of mouse V1 mag. factor
% rho_deg = 0.5;    % (deg / grid)
% smaller value means finer grid, but smaller total range
% too large rho_deg: later W_func_n, K_func_n won't stop decay at max n.
% if you regard 'grid' as 'neuron', then rho_mag * rho_deg should be
  % the spatial density of neuron, which is 10 (um / neuron) now -- fine.


% Adjustable parameters
% alpha_width_C = 1;
% alpha_width_K = 1;
% sigma_A_input_deg = 3.5;    % ~ 6 deg from data
sigma_A_input = sigma_A_input_deg / rho_deg;


%% Input correlation C
DoG_func = @(d, par)...    % (A_e, A_i, sigma_e, sigma_i)    % sigma should be in grid
    par(1) * exp(- d .^ 2 / (2 * par(3) ^ 2)) -...
    par(2) * exp(- d .^ 2 / (2 * par(4) ^ 2));
%
par_DoG_C = [[1.493, 0.495], alpha_width_C * [0.893, 1.587] / rho_deg];
%
C_func = DoG_func(r_grid, par_DoG_C);
%
C_func_n = real(fftshift(fft2(ifftshift(C_func))));


%% Recurrent wiring K
n_r_grid = r_grid / N_sqrt;    % wavenumbers in grid -- n * x must have unit 1;
%
% fft for 2D: see 'supp_How_to_use_fftshift_oddN.m'
% fft2_odd_ctr = @(A) fftshift(fft2(ifftshift(A)));
% ifft2_odd_ctr = @(Z) fftshift(ifft2(ifftshift(Z)));
%
Gn_double = @(n, kappa, sigma_n, sigma_b)...
    kappa * exp(- (2 * pi * sigma_n * n) .^ 2 / 2) +...
    (1 - kappa) * exp(- (2 * pi * sigma_b * n) .^ 2 / 2);
%
sigma_grid = [[5 25] * alpha_width_K; 145, 110]' / (rho_mag * rho_deg);  % row E/I, col N/B
kappa = [0.1, 0.1];
A_rec = [1.05, -1.1] * 5;
% sigma_grid = [[12.5 30] * alpha_width_K; 145, 110]' / (rho_mag * rho_deg);  % row E/I, col N/B
%
W_func_n = A_rec(1) * Gn_double(n_r_grid, kappa(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(n_r_grid, kappa(2), sigma_grid(2, 1), sigma_grid(2, 2));
K_func_n = 1 ./ (1 - W_func_n);
K0_func_n = W_func_n ./ (1 - W_func_n);
K0_func = real(fftshift(ifft2(ifftshift(K0_func_n))));
K_func = K0_func; K_func(idx_ctr, idx_ctr) = K_func(idx_ctr, idx_ctr) + 1;
%
% (Also see in supp_math_FT_iFT_of_recurrent_K.jpg)
% r(x) = (J * r)(x) + h(x) --FT--> r(n) = J(n) * r(n) + h(n)    % h(x) = ∑_β W(x, β)u(β)
  % so, r(n) = h(n) * 1 ./ (1 - J(n)); Define K(n) = 1 ./ (1 - J(n)).
% This requires 2 things:
  % 1) |J(n)| < 1 -- which is fine for our case with E & I
    % (if |J(n)| is small theee are 1 / (1 - J(n)) ~ 1 + J(n), but we don't need that)
  % 2) K(n) is not 'valid' since it asymp. to 1 rather than 0. Cannot integrate.
    % So how to define K(Δx) from Fourier space???
    % Define K0(n) = K(n) - 1 = J(n) ./ (1 - J(n)), which can be iFT to K0(Δx).
    % And K(n) = K0(n) + 1 -> K(Δx) = K0(Δx) + δ(Δx). Consider grid density dx = 1,
      % we can define K(Δx) = {K0(Δx), when Δx ≠ 0} | {K0(Δx) + 1, when Δx = 0}
    % So K(n) and K(Δx) becomes valid when we do FT/iFT with grid density dx = 1.
    % And this is also valid when we consider the problem in convolution form:
      % r(n) = K(n) * h(n) -> r(n) = K0(n) * h(n) + h(n)
      % -- iFT --> r(x) = ∑_y K0(x - y) h(y) + h(x)
      % -> r(x) = ∑_y≠x K0(x - y) h(y) + (K0(0) + 1) h(x)
      % -> r(x) = ∑_y≠x K(x - y) h(y) + K(0) h(x)  % using previous definition of K(Δx)
      % -> r(x) = ∑_y K(x - y) h(y), so now we match the equation r(n) = K(n) * h(n).
    % i.e. the previous way to define K(Δx) is correct!


%% Convert K_func and C_func to translational invariant matrices -- periodic BC of course.
% C_mat = NaN(N, N);
% K_mat = NaN(N, N);
% for k = 1: N
%     i = mod(k - 1, N_sqrt) + 1; j = ceil(k / N_sqrt);
%     ctmp = circshift(C_func, [i - idx_ctr, j - idx_ctr]);
%     C_mat(k, :) = ctmp(:)';
%     ktmp = circshift(K_func, [i - idx_ctr, j - idx_ctr]);
%     K_mat(k, :) = ktmp(:)';
% end
% clear k i j ctmp ktmp
C_mat = circmat_align_2D(repmat(C_func(:)', [N 1]), 0);
K_mat = circmat_align_2D(repmat(K_func(:)', [N 1]), 0);
% tmp = C_mat - DoG_func(get_circ_dist(N_sqrt), par_DoG_C); sum(tmp(:) .^ 2)


%% Arbor function A
% [alpha_j_grid, alpha_i_grid] = meshgrid([-(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2] / N_sqrt);    % normalized to 1
% %
% % target a center on [N_sqrt, N_sqrt] plate -- any point.
% i0 = 1; j0 = 2;
% alpha_j0 = alpha_j_grid(i0, j0); alpha_i0 = alpha_i_grid(i0, j0);
% %
% % get naive dist and circ dist to that target
% dist0 = @(x0, y0) sqrt((alpha_j_grid - x0).^2 + (alpha_i_grid - y0).^2);
% dist_raw = min(cat(3,...
%     dist0(alpha_j0, alpha_i0), dist0(alpha_j0, alpha_i0 - 1), dist0(alpha_j0, alpha_i0 + 1),...
%     dist0(alpha_j0 - 1, alpha_i0), dist0(alpha_j0 - 1, alpha_i0 - 1), dist0(alpha_j0 - 1, alpha_i0 + 1),...
%     dist0(alpha_j0 + 1, alpha_i0), dist0(alpha_j0 + 1, alpha_i0 - 1), dist0(alpha_j0 + 1, alpha_i0 + 1)),...
%     [], 3) * N_sqrt;
% %
% % circshift from that target to everywhere
% dist_circ = NaN(N, N);    % x as row, alpha as column
% for k = 1: N    % each row -- k / (i, j) on 2D x space.
%     i = mod(k - 1, N_sqrt) + 1; j = ceil(k / N_sqrt);
%     dist_tmp = circshift(dist_raw, [i - i0, j - j0]);
%     dist_circ(:, k) = dist_tmp(:);
% end
% clear alpha_j_grid alpha_i_grid i0 j0 dist0 dist_raw k i j dist_tmp
%
dist_circ = circmat_align_2D(repmat(r_grid(:)', [N 1]), 0);
A_mat_in = exp(- dist_circ .^ 2 / (2 * sigma_A_input ^ 2));


%% Initial condition
% Use correct way of L2 norm -- W(:) (vec(W)) is the variable.
L2_norm_func = @(W, L) W * L / sqrt(sum(W(:) .^ 2));
L = 1;
W0 = L2_norm_func(A_mat_in .* (rand([N, N]) - 0.5) * 2, L);
%
eta = 1.5e-4;
% N_step = 8000;
%
% idx_trace = 0: 200: N_step;
N_trace = length(idx_trace);
W_trace = NaN([N, N, N_trace]);
W_trace(:, :, 1) = W0; k_trace = 2;


%% Run
W = W0;
if ifDisplay == 1, fprintf('Start.\n'); tic; end
for ti = 1: N_step
    W = L2_norm_func(W + eta * A_mat_in .* (K_mat * W * C_mat), L);
    if ti == idx_trace(k_trace)
        W_trace(:, :, k_trace) = W;
        k_trace = k_trace + 1;
    end
    if (mod(ti, 200) == 0) & (ifDisplay == 1)
        fprintf([num2str(ti), ' / ', num2str(N_step),...
            ', ', num2str(toc / 60, '%.2f'), ' min.\n']);
    end
end
%
W_trace = reshape(W_trace, [N, N, N_trace]);
clear W k_trace ti

if ifUse73Save == 0
    save(savename, 'N_sqrt', 'N', 'x_sqrt', 'j_grid', 'i_grid', 'r_grid', 'n_r_grid', 'idx_ctr',...
    'rho_mag', 'rho_deg', 'alpha_width_C', 'alpha_width_K', 'sigma_A_input_deg', 'sigma_A_input',...
    'DoG_func', 'par_DoG_C', 'C_func', 'C_func_n', 'C_mat', 'Gn_double', 'sigma_grid', 'kappa', 'A_rec',...
    'W_func_n', 'K_func_n', 'K0_func_n', 'K0_func', 'K_func', 'K_mat', 'dist_circ', 'A_mat_in',...
    'L2_norm_func', 'L', 'eta', 'idx_trace', 'N_step', 'N_trace', 'W_trace');
elseif ifUse73Save == 1
    save(savename, 'N_sqrt', 'N', 'x_sqrt', 'j_grid', 'i_grid', 'r_grid', 'n_r_grid', 'idx_ctr',...
    'rho_mag', 'rho_deg', 'alpha_width_C', 'alpha_width_K', 'sigma_A_input_deg', 'sigma_A_input',...
    'DoG_func', 'par_DoG_C', 'C_func', 'C_func_n', 'C_mat', 'Gn_double', 'sigma_grid', 'kappa', 'A_rec',...
    'W_func_n', 'K_func_n', 'K0_func_n', 'K0_func', 'K_func', 'K_mat', 'dist_circ', 'A_mat_in',...
    'L2_norm_func', 'L', 'eta', 'idx_trace', 'N_step', 'N_trace', 'W_trace', '-v7.3');
end

if ifDisplay == 1, fprintf(['End, ', num2str(toc / 60, '%.2f'), ' min.\n']); end

