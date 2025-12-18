% Skipped numerical eigfunc and fitting / guessing -- in sim_2D_analysis_eig_guess.m
% Power iteration finely tune does not have to be comprehensive,
  % in correct order -- our theoretical solution is accurate enough, if l not too large.

dir_0 = 'D:/Study/CompNeuro/Projects/L4_Oja';
dir_2D = [dir_0, '/sim_theory_2D_fin'];
addpath(dir_0); addpath(dir_2D);
addpath(genpath('D:/Study/CompNeuro/Projects/Functions_simul'));

N_sqrt = 21; N = N_sqrt ^ 2;
rho_deg = 0.5;
alpha_width_C = 1;
alpha_width_K = 1;
sigma_A_input_deg = 3.5;
ifUse73Save = 0;

Ncl = 255; NclH = (Ncl + 1) / 2;
cb = [0 0 0.75]; cm = [1 1 1]; ct = [0.75 0 0];
clr_b2m = [linspace(cb(1), cm(1), NclH)', linspace(cb(2), cm(2), NclH)', linspace(cb(3), cm(3), NclH)'];
clr_m2t = [linspace(cm(1), ct(1), NclH)', linspace(cm(2), ct(2), NclH)', linspace(cm(3), ct(3), NclH)'];
Clr_RWB = [clr_b2m(1: end - 1, :); clr_m2t];
clear Ncl NclH cb cm ct clr_b2m clr_m2t

filename = ['N_sqrt_', num2str(N_sqrt),...
    '_rhodeg', num2str(rho_deg), '_sA', num2str(sigma_A_input_deg), '_aC',...
    num2str(alpha_width_C), '_aK', num2str(alpha_width_K)];



%% Input functions
savename = [dir_2D, '/results/Results_', filename, '.mat'];
%
load(savename, 'par_DoG_C', 'sigma_A_input', 'sigma_grid', 'kappa', 'A_rec', 'Gn_double');
sigma_A = sigma_A_input; 
C1 = par_DoG_C(1); C2 = par_DoG_C(2);
sigma_C_E = par_DoG_C(3); sigma_C_I = par_DoG_C(4);
kappa_K_func = kappa;
clear sigma_A_input par_DoG_C kappa
% A
A0 = 1;
Lambda_A = 1 / (2 * pi * sigma_A);
% C
rho_0 = sqrt( log((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4)) /...
    (2 * (pi^2) * (sigma_C_E ^ 2 - sigma_C_I ^ 2)) );
Cn_max = (2 * pi * C1) *...
    ( (sigma_C_I ^ 2 - sigma_C_E ^ 2) / (sigma_C_I ^ 2 / sigma_C_E ^ 2) )...
    * ( ((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4)) .^...
    (sigma_C_E ^ 2 / (sigma_C_I ^ 2 - sigma_C_E ^ 2)) );
Lambda_C = sqrt( (sigma_C_E ^ 2 - sigma_C_I ^ 2) /...
    (8 * (pi^2) * (sigma_C_E^2) * (sigma_C_I^2) *...
    log((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4))) );
% K -- numerical.
n_fine = (0: 0.0001: N_sqrt) / N_sqrt; % dn = n(2) - n(1);    % (N_sqrt - 1) / 2
Wn = A_rec(1) * Gn_double(n_fine, kappa_K_func(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(n_fine, kappa_K_func(2), sigma_grid(2, 1), sigma_grid(2, 2));
Kn = 1 ./ (1 - Wn);
clear Gn_double A_rec kappa_K_func sigma_grid
% figure; plot(n * N_sqrt, Kn);
[Kn_max, idx] = max(Kn); pK = n_fine(idx); clear idx
%
k_mod_thr = [pK - rho_0, pK + rho_0] * N_sqrt;


%% Other parameters and preparations
% Coordinate
[j_grid, i_grid] = meshgrid((-(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2));
r_grid = sqrt(j_grid .^ 2 + i_grid .^ 2);
theta_grid = atan2(i_grid, j_grid);
idx_ctr = (N_sqrt + 1) / 2; kij_max = idx_ctr - 1;
%
% About R(\rho) solution: Gaussian envelope
Lambda_R = Lambda_A * sqrt( (1 + sqrt( 4 * (Lambda_C^2) / (Lambda_A^2) + 1 )) / 2 );
sigma_R = 1 / (2 * pi * Lambda_R);
Gau_Env_func = exp(- r_grid .^ 2 / (2 * sigma_R ^ 2));
lambda_coeff = (2 * pi * sigma_A^2 * A0) * Cn_max *...
    ( sqrt(2 * pi) * 2 * rho_0 / sqrt(Lambda_A ^ (-2) + Lambda_C ^ (-2) + Lambda_R ^ (-2)) );
%
% Circular Harmonic basis & Gaussian envelope
BslJ_E_cpx = @(m)...
    besselj(abs(m), 2 * pi * rho_0 * r_grid) .* exp(sqrt(-1) * m * theta_grid) .* Gau_Env_func;
m_list = -4: 4; Nm_list = length(m_list);    % 9
BslJ_E_base_cpx = NaN(N, Nm_list);
for m = 1: Nm_list, BslJ_E_base_cpx(:, m) = reshape(BslJ_E_cpx(m_list(m)), [N 1]); end; clear m
%
% Related to solving g_k(\phi) 
dphi_1D = 2 * pi / 360;
phi_1D = 0: dphi_1D: 2 * pi; phi_1D = phi_1D(1: end - 1);
%
tau = Lambda_A / rho_0;
a_func = @(delta_phi) exp(- delta_phi .^ 2 / (2 * tau ^ 2));
f_wrap_phi = @(d_phi) atan2(sin(d_phi), cos(d_phi));
a_mat = a_func(f_wrap_phi(bsxfun(@minus, phi_1D, phi_1D')));
%
exp_neg_i_m_phi = exp(-sqrt(-1) * m_list' * phi_1D);
%
N_deg_half = 7;    % i.e. l = 0, 1, 2, ..., 6
N_deg_full = N_deg_half * 2;    % consider conjugated pair of vec{k} and -\vec{k} => Take Re and Im.
%
% Related to finely tune
load(savename, 'A_mat_in', 'K_mat', 'C_mat');
N_itr_max = 150;
eigval_itr_thr = 1e-5;
L2_norm = @(W) W / sqrt(sum(W(:) .^ 2));



%% HO theory of eigval
k_mod_theory = (0: 0.001: N_sqrt) / N_sqrt;
%
k_mod_type = NaN(size(k_mod_theory));
k_mod_type(k_mod_theory < pK - rho_0) = 1;
k_mod_type((k_mod_theory >= pK - rho_0) & (k_mod_theory <= pK + rho_0)) = 2;
k_mod_type(k_mod_theory > pK + rho_0) = 3;
%
f_star = NaN(size(k_mod_theory));
f_star(k_mod_type == 1) = interp1(n_fine, Kn, rho_0 + k_mod_theory(k_mod_type == 1));
f_star(k_mod_type == 2) = Kn_max;
f_star(k_mod_type == 3) = interp1(n_fine, Kn, k_mod_theory(k_mod_type == 3) - rho_0);
%
Kn_p1 = gradient(Kn, n_fine);
Kn_p2 = gradient(Kn_p1, n_fine);
%
kappa = NaN(size(k_mod_theory));
kappa(k_mod_type == 1) = interp1(n_fine, Kn_p1, rho_0 + k_mod_theory(k_mod_type == 1))...
    .* ( (rho_0 * k_mod_theory(k_mod_type == 1)) ./ (rho_0 + k_mod_theory(k_mod_type == 1)) );
kappa(k_mod_type == 2) = ( interp1(n_fine, Kn_p2, pK) / (4 * (pK^2)) )...
    * ((rho_0 + k_mod_theory(k_mod_type == 2)) .^ 2 - pK ^ 2)...
    .* ((rho_0 - k_mod_theory(k_mod_type == 2)) .^ 2 - pK ^ 2);
kappa(k_mod_type == 3) = -interp1(n_fine, Kn_p1, abs(rho_0 - k_mod_theory(k_mod_type == 3)))...
    .* ( (rho_0 * k_mod_theory(k_mod_type == 3)) ./ (k_mod_theory(k_mod_type == 3) - rho_0) );
%
lambda_0 = (4 * pi^2 * A0 * Cn_max) /...
    sqrt(1 / (Lambda_A^2) + 1 / (Lambda_R^2));    % + 1 / (Lambda_C^2)
Eigval_theory = lambda_0 * (f_star - (tau / 2) * sqrt(f_star .* kappa));    %  + (tau^2)/4 * kappa
% eigenvalue0(:, 1) = lambda_0 * interp1(n_fine, Kn, rho_0);
%
% figure; plot(k_mod_theory * N_sqrt, Eigval_theory);


%% Numerical accurate eigs -- numerical g_k(\phi) & finely tune
% How to solve angular g_k(\phi).
% \mu = \lambda / (\widetilde{C}(\rho_0) rho_0 \sqrt(2 * pi) Lambda_RC)
% \mu g_k(\phi) = \int_0^{2 \pi} a(\phi - \phi') f_k(\phi') g_k(\phi') d\phi    (*)
% where
% a(\Delta\phi) = (2 \pi \sigma_A^2 A0) exp(-\Delta\phi^2 / (2 \tau^2))
% f_k(\phi) = \widetilde{K}(\sqrt{\rho_0^2 + |k|^2 - 2 \rho_0 |k| cos(\phi - \phi_k)})
%
% Convert (*) to matrix -- \mu * g_k = ((a_ij f_j) * d\phi) * g_k
% Let g_k = g0_k ./ sqrt(f) -- \mu * g0_k = ((sqrt(f)_i a_ij sqrt(f)_j) * d\phi) * g0_k
%
% This numerical eigval and lambda_factor gets wrong -- almost double of final real eigval.
% But anyway -- accurate eigval need finely tune &
  % quali peak mode conclusion not affected. So skip it here.
% Here we just use numerical g_k(\phi) for finely tune ic
  % -- they're accurate enough so no need re-order, if l not too large.

R_theory = cell(1, (N + 1) / 2);
coeff_theory = cell(1, (N + 1) / 2);
R_theory_FinelyTuned = cell(1, (N + 1) / 2);
Eigval_theory_FinelyTuned = cell(1, (N + 1) / 2);
Eigval_theory_FinelyTuned_Itr = cell(1, (N + 1) / 2);
R_theory_FinelyTuned_Rsqr_Itr = cell(1, (N + 1) / 2);
Eig_if_valid = cell(1, (N + 1) / 2);
%
tic;
kji_idx = 0;
for kj = 0: kij_max
if kj == 0, ki_start = 0; else, ki_start = -kij_max; end
for ki = ki_start: kij_max
kji_idx = kji_idx + 1;
% kj = 1; ki = 1;
%
k_mod = sqrt(sum([kj, ki] .^ 2)) / N_sqrt;
phi_k = atan2(ki, kj);
exp_kx = exp(2 * pi * sqrt(-1) * (kj * j_grid(:) + ki * i_grid(:)) / N_sqrt);
%
% Numerical solution of g_k(\phi), and "theoretical" eigvec for fine tune ic.
if ~(k_mod == 0)
    f_k_phi = interp1(n_fine, Kn,...
        sqrt(rho_0 ^ 2 + k_mod ^ 2 - 2 * rho_0 * k_mod * cos(phi_1D - phi_k)));
    % figure; hold on;
    % plot(phi_1D, f_k_phi);
    % xlim([0 2 * pi]); set(gca, 'XTick', linspace(0, 2 * pi, 9));
    %
    % plot(pi * ones(1, 2), [1.1 1.8], 'k--');
    % [~, idxM] = max(f_k_phi); plot(phi_1D(idxM) * ones(1, 2), [1.1 1.8], 'r--');
    % plot(2 * pi - phi_1D(idxM) * ones(1, 2), [1.1 1.8], 'r--');
    % exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/f_k_phi_3_0.eps'],...
    %     'BackgroundColor', 'none', 'ContentType', 'vector');
    %
    Operator_self_adjoint = (diag(sqrt(f_k_phi)) * a_mat * diag(sqrt(f_k_phi))) * dphi_1D;
    [U_g0_k, S_mu, ~] = svd(Operator_self_adjoint);
    gk_phi_1D = bsxfun(@times, 1 ./ sqrt(f_k_phi'), U_g0_k(:, 1: N_deg_half));
    % S_mu = diag(S_mu); lambda_theo = S_mu(1: N_deg_half) * lambda_coeff;
    clear Operator_self_adjoint U_g0_k S_mu
    % [V, D] = eig(a_mat * diag(f_k_phi) * dphi_1D); diag(D)
    %
    Rk_theory = NaN(N ^ 2, N_deg_full);
    ckl_theory = NaN(Nm_list, N_deg_full);
    for l = 1: N_deg_half
        ckl_hat = (sqrt(-1) .^ abs(m_list')) .*...
            sum(bsxfun(@times, gk_phi_1D(:, l).', exp_neg_i_m_phi) * dphi_1D, 2);    % (9, 1)
        % Wkl_hat = BslJ_E_base_cpx * ckl_hat;
        % figure;
        % subplot(1, 2, 1); imagesc(reshape(real(Wkl_hat), [N_sqrt N_sqrt])); axis square; colorbar;
        % subplot(1, 2, 2); imagesc(reshape(imag(Wkl_hat), [N_sqrt N_sqrt])); axis square; colorbar;
        Wkl_cpx = exp_kx * (BslJ_E_base_cpx * ckl_hat).';
        ckl_theory(:, [2 * l - 1, 2 * l]) = repmat(ckl_hat, [1 2]);
        Rk_theory(:, 2 * l - 1) =...
            reshape(circmat_align_2D(L2_norm(real(Wkl_cpx)), 0), [N^2 1]);
        Rk_theory(:, 2 * l) =...
            reshape(circmat_align_2D(L2_norm(imag(Wkl_cpx)), 0), [N^2 1]);
    end
    % Check shape of f_k_phi, gk_phi_1D!
    clear f_k_phi gk_phi_1D l ckl_hat Wkl_cpx
elseif k_mod == 0
    gk_phi_1D = NaN(length(phi_1D), N_deg_half);
    for l = 0: (N_deg_half - 1)
        gk_phi_1D(:, l + 1) = exp(sqrt(-1) * l * phi_1D);
    end
    %
    Rk_theory = NaN(N ^ 2, N_deg_half);
    ckl_theory = NaN(Nm_list, N_deg_half);
    for l = 1: N_deg_half
        ckl_hat = (sqrt(-1) .^ abs(m_list')) .*...
            sum(bsxfun(@times, gk_phi_1D(:, l).', exp_neg_i_m_phi) * dphi_1D, 2);    % (9, 1)
        Wkl_cpx = exp_kx * (BslJ_E_base_cpx * ckl_hat).';
        ckl_theory(:, l) = ckl_hat;
        Rk_theory(:, l) =...
            reshape(circmat_align_2D(L2_norm(real(Wkl_cpx)), 0), [N^2 1]);
    end
    % Check shape of f_k_phi, gk_phi_1D!
    clear f_k_phi gk_phi_1D l ckl_hat Wkl_cpx
end
%
% Finely Tune
Qk_theory_FinelyTuned = NaN(size(Rk_theory));
Rk_theory_FinelyTuned = NaN(size(Rk_theory));
Eigval_k_theory_FinelyTuned = NaN(1, size(Rk_theory, 2));
Eigval_k_theory_FinelyTuned_Itr = NaN(N_itr_max + 1, size(Rk_theory, 2));
Rk_theory_FinelyTuned_Rsqr_Itr = NaN(N_itr_max + 1, size(Rk_theory, 2));
for l = 1: size(Rk_theory, 2)
    [Qk_theory_FinelyTuned(:, l), Rk_theory_FinelyTuned(:, l),...
        Eigval_k_theory_FinelyTuned(l), ~, Eigval_k_theory_FinelyTuned_Itr(:, l),...
        Rk_theory_FinelyTuned_Rsqr_Itr(:, l)] =...
    func_eigvec_FinelyTune(Rk_theory(:, l),...
        Qk_theory_FinelyTuned(:, 1: l - 1),...
        A_mat_in, K_mat, C_mat, N_itr_max, eigval_itr_thr);
    % Just use previous l of same k as Q_history
end
clear Qk_theory_FinelyTuned l
%
% l = 1;
% RF_2D_visualize_fast(circmat_align_2D(reshape(Rk_theory(:, l),...
%     [N N]), 1), 1, -0.015, 0.015, Clr_RWB);
% RF_2D_visualize_fast(circmat_align_2D(reshape(Rk_theory_FinelyTuned(:, l),...
%     [N N]), 1), 1, -0.015, 0.015, Clr_RWB);
% figure;
% subplot(2, 1, 1); plot(Eigval_k_theory_FinelyTuned_Itr(:, l));
% subplot(2, 1, 2); plot(Rk_theory_FinelyTuned_Rsqr_Itr(:, l));
%
Eig_if_valid_k = ones(size(Rk_theory, 2), 1);
for l = 1: size(Rk_theory, 2)
    x = 0: N_itr_max;
    d2y = gradient(gradient(Eigval_k_theory_FinelyTuned_Itr(:, l)', x), x);
    if any(d2y(5: end - 1) .* d2y(6: end) < 0)
        Eig_if_valid_k(l: end) = 0; break;
    end
end
clear jl x d2y
%
R_theory{kji_idx} = Rk_theory;
coeff_theory{kji_idx} = ckl_theory;
R_theory_FinelyTuned{kji_idx} = Rk_theory_FinelyTuned;
Eigval_theory_FinelyTuned{kji_idx} = Eigval_k_theory_FinelyTuned;
Eigval_theory_FinelyTuned_Itr{kji_idx} = Eigval_k_theory_FinelyTuned_Itr;
R_theory_FinelyTuned_Rsqr_Itr{kji_idx} = Rk_theory_FinelyTuned_Rsqr_Itr;
Eig_if_valid{kji_idx} = Eig_if_valid_k;
%
fprintf(['k (', num2str(kj), ', ', num2str(ki), '), ', num2str(kji_idx), ' / ',...
    num2str((N + 1) / 2), ' done,  Finely tune ', num2str(sum(Eig_if_valid_k)),...
    ' / ', num2str(size(Rk_theory, 2)), ' good, ', num2str(toc / 60, '%.2f'), ' min.\n']);
end
end
%
clear kj ki kji_idx k_mod k_phi exp_kx
clear Rk_theory ckl_theory Rk_theory_FinelyTuned Eigval_k_theory_FinelyTuned
clear Eigval_k_theory_FinelyTuned_Itr Rk_theory_FinelyTuned_Rsqr_Itr Eig_if_valid_k
%
% save([dir_2D, '/results/Eigs_HalfTheory_FinelyTuned_', filename, '.mat'],...
%     'R_theory', 'coeff_theory', 'R_theory_FinelyTuned', 'Eigval_theory_FinelyTuned',...
%     'Eigval_theory_FinelyTuned_Itr', 'R_theory_FinelyTuned_Rsqr_Itr',...
%     'Eig_if_valid', 'kij_max', '-v7.3');

% Reorganize
N_eigvec_valid = 0;
for kji_idx = 1: (N + 1) / 2
    N_eigvec_valid = N_eigvec_valid + sum(Eig_if_valid{kji_idx});
end
%
level_valid = NaN(N_eigvec_valid, 1);
R_theory_valid = NaN(N ^ 2, N_eigvec_valid);
coeff_theory_valid = NaN(Nm_list, N_eigvec_valid);
R_FinelyTuned_valid = NaN(N ^ 2, N_eigvec_valid);
kj_ki_valid = NaN(N_eigvec_valid, 2);
Eigval_FinelyTuned_valid = NaN(N_eigvec_valid, 1);
Eigval_FinelyTuned_Itr_valid = NaN(N_itr_max + 1, N_eigvec_valid);
Rsqr_FinelyTuned_Itr_valid = NaN(N_itr_max + 1, N_eigvec_valid);
%
kji_idx = 0; k_valid = 0;
for kj = 0: kij_max
    if kj == 0, ki_start = 0; else, ki_start = -kij_max; end
    for ki = ki_start: kij_max
        kji_idx = kji_idx + 1;
        %
        idx_l = find(Eig_if_valid{kji_idx} == 1);
        idx = k_valid + (1: length(idx_l));
        %
        if (kj == 0) & (ki == 0), level_valid(idx) = idx_l - 1;
        else, level_valid(idx) = ceil(idx_l / 2) - 1; end
        kj_ki_valid(idx, :) = repmat([kj, ki], [length(idx_l) 1]);
        R_theory_valid(:, idx) = R_theory{kji_idx}(:, idx_l);
        coeff_theory_valid(:, idx) = coeff_theory{kji_idx}(:, idx_l);
        R_FinelyTuned_valid(:, idx) = R_theory_FinelyTuned{kji_idx}(:, idx_l);
        Eigval_FinelyTuned_valid(idx) = Eigval_theory_FinelyTuned{kji_idx}(idx_l);
        Eigval_FinelyTuned_Itr_valid(:, idx) = Eigval_theory_FinelyTuned_Itr{kji_idx}(:, idx_l);
        Rsqr_FinelyTuned_Itr_valid(:, idx) = R_theory_FinelyTuned_Rsqr_Itr{kji_idx}(:, idx_l);
        %
        k_valid = k_valid + length(idx_l);
        %
        % save space
        R_theory{kji_idx} = [];
        R_theory_FinelyTuned{kji_idx} = [];
    end
end
clear R_theory coeff_theory R_theory_FinelyTuned Eigval_theory_FinelyTuned
clear Eigval_theory_FinelyTuned_Itr R_theory_FinelyTuned_Rsqr_Itr
%
[Eigval_FinelyTuned_valid, order] = sort(Eigval_FinelyTuned_valid, 'descend');
level_valid = level_valid(order);
kj_ki_valid = kj_ki_valid(order, :);
R_theory_valid = R_theory_valid(:, order);
coeff_theory_valid = coeff_theory_valid(:, order);
R_FinelyTuned_valid = R_FinelyTuned_valid(:, order);
Eigval_FinelyTuned_Itr_valid = Eigval_FinelyTuned_Itr_valid(:, order);
Rsqr_FinelyTuned_Itr_valid = Rsqr_FinelyTuned_Itr_valid(:, order);
%
save([dir_2D, '/results/Eigs_HalfTheory_FinelyTuned_', filename, '.mat'],...
    'k_mod_theory', 'Eigval_theory', 'level_valid', 'kj_ki_valid',...
    'R_theory_valid', 'coeff_theory_valid', 'R_FinelyTuned_valid',...
    'Eigval_FinelyTuned_valid', 'Eigval_FinelyTuned_Itr_valid',...
    'Rsqr_FinelyTuned_Itr_valid', 'Eig_if_valid', 'kij_max', '-v7.3');



%% Eigvec plots (see below)
% Plot iteration process of finely tune
% Find out J basis coeff of R for plots.
% Figure / video: coeff & iteration process of each eigvec

% R-sqr & coeff over idx
load([dir_2D, '/results/Eigs_HalfTheory_FinelyTuned_', filename, '.mat'],...
    'kj_ki_valid', 'Eigval_FinelyTuned_valid', 'Rsqr_FinelyTuned_Itr_valid',...
    'R_FinelyTuned_valid', 'R_theory_valid');
N_eig_count = 128;
R_FinelyTuned_valid = R_FinelyTuned_valid(:, 1: N_eig_count);
R_theory_valid = R_theory_valid(:, 1: N_eig_count);

r_grid = sqrt(j_grid .^ 2 + i_grid .^ 2);
theta_grid = atan2(i_grid, j_grid);
Env_Gau_func = exp(- r_grid .^ 2 / (2 * sigma_R ^ 2));
%
n0_r = 2 * pi * rho_0 * r_grid;
m_list = -4: 4; Nm_list = length(m_list);
%
BslJ_E_cos = @(m) besselj(m, n0_r) .* cos(m * theta_grid) .* Env_Gau_func;
BslJ_E_sin = @(m) besselj(m, n0_r) .* sin(m * theta_grid) .* Env_Gau_func;
BslJ_E_base = reshape(cat(3, BslJ_E_cos(0), BslJ_E_cos(1), BslJ_E_sin(1),...
    BslJ_E_cos(2), BslJ_E_sin(2), BslJ_E_cos(3), BslJ_E_sin(3),...
    BslJ_E_cos(4), BslJ_E_sin(4)), [N Nm_list]);
%
cos_kx = @(kj, ki) cos(2 * pi * (kj * j_grid(:) + ki * i_grid(:)) / N_sqrt);
sin_kx = @(kj, ki) sin(2 * pi * (kj * j_grid(:) + ki * i_grid(:)) / N_sqrt);
%
R_sqr_func_v = @(v, v0) 1 - sum((v - v0) .^ 2) / sum((v0 - mean(v0)) .^ 2);


R_sqr_FinelyTuned = NaN(1, N_eig_count);
Coeff_theory = NaN(2 * Nm_list, N_eig_count);
%
for idx_eigvec = 1: N_eig_count
    kj = kj_ki_valid(idx_eigvec, 1); ki = kj_ki_valid(idx_eigvec, 2);
    G_F_B_base = NaN(N ^ 2, 2 * Nm_list);
    % No (0, 0) for first 128.
    for j = 1: Nm_list
        G_F_B_base(:, 2 * j - 1) =...
            L2_norm(reshape(cos_kx(kj, ki) * BslJ_E_base(:, j)', [N^2, 1]));
        G_F_B_base(:, 2 * j) =...
            L2_norm(reshape(sin_kx(kj, ki) * BslJ_E_base(:, j)', [N^2, 1]));
    end
    Rk_ctr_theory = circmat_align_2D(reshape(R_theory_valid(:, idx_eigvec), [N N]), 1);
    Rk_ctr_FinelyTuned = circmat_align_2D(reshape(R_FinelyTuned_valid(:, idx_eigvec), [N N]), 1);
    %
    % % tmp = G_F_B_base * coeff_theory;
    % % R_sqr_theo = R_sqr_func_v(tmp(:), Rk_ctr_theory(:));    % must be 1
    % coeff_FinelyTuned = pinv(G_F_B_base) * Rk_ctr_FinelyTuned(:);
    % tmp1 = G_F_B_base * coeff_theory; tmp2 = G_F_B_base * coeff_FinelyTuned;
    % R2_1(idx_eigvec) = R_sqr_func_v(tmp1(:), Rk_ctr_FinelyTuned(:));
    % R2_2(idx_eigvec) = R_sqr_func_v(tmp2(:), Rk_ctr_FinelyTuned(:));
    %
    Coeff_theory(:, idx_eigvec) = pinv(G_F_B_base) * Rk_ctr_theory(:);
    % Coeff_FinelyTuned(:, idx_eigvec) = pinv(G_F_B_base) * Rk_ctr_FinelyTuned(:);
    Rk_ctr_FinelyTuned_within = G_F_B_base * Coeff_theory(:, idx_eigvec);
    R_sqr_FinelyTuned(idx_eigvec) =...
        R_sqr_func_v(Rk_ctr_FinelyTuned_within(:), Rk_ctr_FinelyTuned(:));
end
clear kj ki G_F_B_base j Rk_ctr_theory Rk_ctr_FinelyTuned Rk_ctr_FinelyTuned_within


figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 0.67]);
subplot(4, 1, 1);
plot(Eigval_FinelyTuned_valid(1: N_eig_count));
ylim([15.5 16.5]);
xlim([1 128]); set(gca, 'XTick', [1 8: 8: 128]);
%
subplot(4, 1, 2); hold on;
plot(R_sqr_FinelyTuned, 'r');
%plot(Rsqr_FinelyTuned_Itr_valid(end, 1: N_eig_count), 'b--');
ylim([0.94 1]);
xlim([1 128]); set(gca, 'XTick', [1 8: 8: 128]);
%
subplot(4, 1, 3);
imagesc(Coeff_theory(1: 2: end, :)); colormap(Clr_RWB); colorbar; clim([-1 1]);
set(gca, 'XTick', [1 8: 8: 128], 'YTick', 1: 9);
%
subplot(4, 1, 4);
imagesc(Coeff_theory(2: 2: end, :)); colormap(Clr_RWB); colorbar; clim([-1 1]);
set(gca, 'XTick', [1 8: 8: 128], 'YTick', 1: 9);
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename,...
    '/eigs/Eigvec_coeff.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/Eigvec_coeff.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;




%% Eigval - |k| !!!
load([dir_2D, '/results/Eigs_HalfTheory_FinelyTuned_', filename, '.mat'],...
    'level_valid', 'kj_ki_valid', 'Eigval_FinelyTuned_valid', 'kij_max',...
    'k_mod_theory', 'Eigval_theory');
k_mod = sqrt(sum(kj_ki_valid .^ 2, 2));
%
l_tot = 6;%length(unique(level_valid));
%clr = jet(l_tot + 1); clr = clr([1: 5, 7], :); clr = clr(end: -1: 1, :) * 0.9;
clr = [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 0.1 0 0.5] * 0.9;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0, 0.5, 0.8]); hold on;
lg = zeros(1, l_tot + 1); lgdtxt = cell(1, l_tot + 1);
lg(1) = plot(k_mod_theory * N_sqrt, Eigval_theory, 'color', clr(1, :), 'linewidth', 3);
lgdtxt{1} = 'l = 0 (theory)';
for l = 0: l_tot - 1
    idx_l = find(level_valid == l);
    lg(l + 2) = scatter(k_mod(idx_l), Eigval_FinelyTuned_valid(idx_l),...
        50, clr(l + 1, :), 'filled');
    lgdtxt{l + 2} = ['l = ', num2str(l)];
end
plot((pK - rho_0) * N_sqrt * ones(1, 2), [9 17], 'k--');
plot(pK * N_sqrt * ones(1, 2), [9 17], 'k--');
plot((pK + rho_0) * N_sqrt * ones(1, 2), [9 17], 'k--');
text((pK - rho_0) * N_sqrt, 17.2, '$\rho_K^0 - \rho_C^0$', 'fontsize', 12,...
    'HorizontalAlignment', 'center', 'interpreter', 'latex');
text(pK * N_sqrt, 17.2, '$\rho_K^0$', 'fontsize', 12,...
    'HorizontalAlignment', 'center', 'interpreter', 'latex');
text((pK + rho_0) * N_sqrt, 17.2, '$\rho_K^0 + \rho_C^0$', 'fontsize', 12,...
    'HorizontalAlignment', 'center', 'interpreter', 'latex');
legend(lg, lgdtxt);
axis square; axis([0 15 9 17]); set(gca, 'XTick', 0: 1: 14);
xlabel('$|\vec{k}|$', 'interpreter', 'latex');
ylabel('Eigenvalue');
ax = gca; ax.FontSize = 16;
%
pause(0.5);
print(gcf, '-dpng', [dir_2D, '/figures/', filename,...
    '/eigs/Eigval_2D.png']);


load(savename, 'rho_mag', 'rho_deg');
d_unit = (rho_mag * rho_deg);
xmax = 50;    % 1000 / 21
%
l_tot = 6;%length(unique(level_valid));
%clr = jet(l_tot + 1); clr = clr([1: 5, 7], :); clr = clr(end: -1: 1, :) * 0.9;
clr = [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 0.1 0 0.5] * 0.9;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0, 0.5, 0.8]); hold on;
clear lg lgdtxt
lg = zeros(1, l_tot + 1); lgdtxt = cell(1, l_tot + 1);
lg(1) = plot(k_mod_theory * (1000 / d_unit), Eigval_theory, 'color', clr(1, :), 'linewidth', 3);
lgdtxt{1} = 'l = 0 (theory)';
for l = 0: l_tot - 1
    idx_l = find(level_valid == l);
    lg(l + 2) = scatter((k_mod(idx_l) / N_sqrt) * (1000 / d_unit),...
        Eigval_FinelyTuned_valid(idx_l), 50, clr(l + 1, :), 'filled');
    lgdtxt{l + 2} = ['l = ', num2str(l)];
end
plot((pK - rho_0) * (1000 / d_unit) * ones(1, 2), [9 17], 'k--');
plot((pK + rho_0) * (1000 / d_unit) * ones(1, 2), [9 17], 'k--');
text((pK - rho_0) * (1000 / d_unit), 17.2, '$\rho_K^0 - \rho_C^0$', 'fontsize', 12,...
    'HorizontalAlignment', 'center', 'interpreter', 'latex');
text((pK + rho_0) * (1000 / d_unit), 17.2, '$\rho_K^0 + \rho_C^0$', 'fontsize', 12,...
    'HorizontalAlignment', 'center', 'interpreter', 'latex');
legend(lg, lgdtxt);
axis square; axis([0 xmax 9 17]); set(gca, 'XTick', 0: 5: xmax);
xlabel('$|\vec{k}|$', 'interpreter', 'latex');
ylabel('Eigenvalue');
ax = gca; ax.FontSize = 16;
%
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/Eigval_2D.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;






%% Schematic: Plot circular harmonic functions
r_grid = sqrt(j_grid .^ 2 + i_grid .^ 2);
theta_grid = atan2(i_grid, j_grid);
Env_Gau_func = exp(- r_grid .^ 2 / (2 * sigma_R ^ 2));
%
n0_r = 2 * pi * rho_0 * r_grid;
% BslJ_E_cpx = @(m) besselj(abs(m), n0_r) .* exp(sqrt(-1) * m * theta_grid) .* Env_Gau_func;
m_list = -4: 4; Nm_list = length(m_list);
%
BslJ_E_cos = @(m) besselj(m, n0_r) .* cos(m * theta_grid) .* Env_Gau_func;
BslJ_E_sin = @(m) besselj(m, n0_r) .* sin(m * theta_grid) .* Env_Gau_func;
BslJ_E_base = reshape(cat(3, BslJ_E_cos(0), BslJ_E_cos(1), BslJ_E_sin(1),...
    BslJ_E_cos(2), BslJ_E_sin(2), BslJ_E_cos(3), BslJ_E_sin(3),...
    BslJ_E_cos(4), BslJ_E_sin(4)), [N Nm_list]);
%
BslJ_cos = @(m) besselj(m, n0_r) .* cos(m * theta_grid);
BslJ_sin = @(m) besselj(m, n0_r) .* sin(m * theta_grid);
BslJ_base = reshape(cat(3, BslJ_cos(0), BslJ_cos(1), BslJ_sin(1),...
    BslJ_cos(2), BslJ_sin(2), BslJ_cos(3), BslJ_sin(3),...
    BslJ_cos(4), BslJ_sin(4)), [N Nm_list]);
%
G_F_B_title_txt = {'0(2\pi\rho_0r)',...
    '1(2\pi\rho_0r) \cos(\theta)', '1(2\pi\rho_0r) \sin(\theta)',...
    '2(2\pi\rho_0r) \cos(2\theta)', '2(2\pi\rho_0r) \sin(2\theta)',...
    '3(2\pi\rho_0r) \cos(3\theta)', '3(2\pi\rho_0r) \sin(3\theta)',...
    '4(2\pi\rho_0r) \cos(4\theta)', '4(2\pi\rho_0r) \sin(4\theta)'};
G_F_B_clim = [-0.5 1; -0.4 0.4; -0.4 0.4; -0.4 0.4; -0.4 0.4;...
    -0.2 0.2; -0.2 0.2; -0.1 0.1; -0.1 0.1];
xtk = 1: 2: N_sqrt; xtk_lb = xtk - ((N_sqrt + 1) / 2);
plot_idx = [1: 5, 7: 10];
%
for i = 1: 2
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05, 0, 1, 0.6]);
for j = 1: length(plot_idx)
    plot_k = plot_idx(j);
    subplot(2, 5, plot_k);
    if i == 1, BB = BslJ_E_base; elseif j == 2, BB = BslJ_base; end
    imagesc(reshape(BB(:, j), [N_sqrt N_sqrt]));
    axis square; colormap(Clr_RWB); colorbar;
    if i == 1, clim(G_F_B_clim(j, :)); elseif j == 2, clim([-0.5 0.5]); end
    set(gca, 'XTick', xtk, 'XTickLabel', xtk_lb, 'YTick', xtk, 'YTickLabel', xtk_lb);
    xlabel('$a_j - x_j$', 'interpreter', 'latex');
    ylabel('$a_i - x_i$', 'interpreter', 'latex');
    if i == 1
        title(['$E(r)J_', G_F_B_title_txt{j}, '$'],...
            'interpreter', 'latex', 'fontweight', 'normal');
        name_tail = '_withE';
    elseif i == 2
        title(['$J_', G_F_B_title_txt{j}, '$'],...
            'interpreter', 'latex', 'fontweight', 'normal');
        name_tail = '_withoutE';
    end
end
pause(0.5);
print(gcf, '-dpng', [dir_2D, '/figures/', filename,...
    '/eigs/Circular_harmonic_bases', name_tail, '.png']);
close;
end

