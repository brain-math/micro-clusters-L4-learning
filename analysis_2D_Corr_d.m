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
savename = [dir_2D, '/results/Results_', filename, '.mat'];



% Load W(t)
load(savename, 'W_trace', 'idx_trace', 'N_trace', 'N_step', 'rho_mag', 'rho_deg', 'r_grid');
d_unit = (rho_mag * rho_deg); clear rho_mag rho_deg    % 10 um
Wt = W_trace(:, :, end);
Wt_ctr = circmat_align_2D(Wt, 1);

% Get input pars -- rho_C, rho_K, sigma_R
load(savename, 'par_DoG_C', 'sigma_A_input', 'sigma_grid', 'kappa', 'A_rec', 'Gn_double');
sigma_A = sigma_A_input; 
C1 = par_DoG_C(1); C2 = par_DoG_C(2);
sigma_C_E = par_DoG_C(3); sigma_C_I = par_DoG_C(4);
kappa_K_func = kappa;
clear sigma_A_input par_DoG_C kappa
%
Lambda_A = 1 / (2 * pi * sigma_A);
rho_C = sqrt( log((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4)) /...
    (2 * (pi^2) * (sigma_C_E ^ 2 - sigma_C_I ^ 2)) );
Lambda_C = sqrt( (sigma_C_E ^ 2 - sigma_C_I ^ 2) /...
    (8 * (pi^2) * (sigma_C_E^2) * (sigma_C_I^2) *...
    log((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4))) );
Lambda_R = Lambda_A * sqrt( (1 + sqrt( 4 * (Lambda_C^2) / (Lambda_A^2) + 1 )) / 2 );
sigma_R = 1 / (2 * pi * Lambda_R);
clear sigma_A C1 C2 sigma_C_E sigma_C_I Lambda_A Lambda_C Lambda_R
%
n_fine = (0: 0.0001: N_sqrt) / N_sqrt; % dn = n(2) - n(1);    % (N_sqrt - 1) / 2
Wn = A_rec(1) * Gn_double(n_fine, kappa_K_func(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(n_fine, kappa_K_func(2), sigma_grid(2, 1), sigma_grid(2, 2));
Kn = 1 ./ (1 - Wn);
[Kn_max, idx] = max(Kn);
rho_K = n_fine(idx);
clear n_fine Wn Gn_double A_rec kappa_K_func sigma_grid idx Kn Kn_max


%% Sim inner product
% % Way 1
% load(savename, 'x_sqrt');
% [j_grid, i_grid] = meshgrid(x_sqrt);
% %
% dtheta_deg = 15;
% theta_stim = [dtheta_deg: dtheta_deg: 180] * (pi / 180); N_theta = length(theta_stim);
% dphi_deg = 15;
% phi_stim = [dphi_deg: dphi_deg: 360] * (pi / 180); N_phi = length(phi_stim);
% n0_list = [0: 0.5: 10] / N_sqrt;
% N_n0 = length(n0_list);
% %
% stim_grating = NaN(N, N_phi, N_theta, N_n0);
% for n0_j = 1: N_n0
%     n0 = n0_list(n0_j);
%     for theta_i = 1: N_theta
%         theta = theta_stim(theta_i);
%         for phi_k = 1: N_phi 
%             stim_grating(:, phi_k, theta_i, n0_j) = cos(2 * pi * n0 *...
%                 (j_grid(:) * cos(theta) + i_grid(:) * sin(theta)) + phi_stim(phi_k));
%         end
%     end
% end
% r_ffwd_sqr = NaN(N, N_theta, N_n0);
% for n0_j = 1: N_n0
% for theta_i = 1: N_theta
%     stim_grating_ij = stim_grating(:, :, theta_i, n0_j);    % (N, N_phi);
%     rL = Wt_ctr * stim_grating_ij;
%     rLp_sqr = (rL .* (rL > 0)) .^ 2;
%     r_ffwd_sqr(:, theta_i, n0_j) = max(rLp_sqr, [], 2);
% end
% end
% clear n0_j n0 theta_i theta phi_k stim_grating_ij rL rLp_sqr j_grid i_grid x_sqrt stim_grating
% %
% r_ffwd_sqr = reshape(r_ffwd_sqr, [N, N_theta * N_n0]);
% r_inner_product_1 = triu_new(r_ffwd_sqr * r_ffwd_sqr', 0, 1);
% clear r_ffwd_sqr

% Way 2
Wt_ctr_absFTsqr_trace = NaN(N, N);
for xi = 1: N
    Wx = reshape(Wt_ctr(xi, :), [N_sqrt, N_sqrt]);
    ft_tmp = abs(fftshift(fft2(ifftshift(Wx))));
    Wt_ctr_absFTsqr_trace(xi, :) = (ft_tmp(:) .^ 2)';
end
clear Wx ft_tmp xi
r_inner_product = triu_new(Wt_ctr_absFTsqr_trace * Wt_ctr_absFTsqr_trace', 0, 1);


% Bin over d
load(savename, 'dist_circ');
d_BinCenter = [1, sqrt(2), 2, sqrt(5), sqrt(8), 3.5: 0.5: 10];
N_BinCenter = length(d_BinCenter);
d_BinEdge = [d_BinCenter(1) - (d_BinCenter(2) - d_BinCenter(1)) / 2,...
    mean([d_BinCenter(1: N_BinCenter - 1); d_BinCenter(2: N_BinCenter)], 1),...
    d_BinCenter(end) + (d_BinCenter(end) - d_BinCenter(end - 1)) / 2];
%
dist_circ = triu_new(dist_circ, 0, 1);
idx_valid = find(dist_circ < d_BinEdge(end));
%
% [r_inner_product_1_avg, r_inner_product_1_se, ~] =...
%     histogram_mean_sem(r_inner_product_1(idx_valid), dist_circ(idx_valid), d_BinEdge);
[r_inner_product_avg, r_inner_product_se, ~] =...
    histogram_mean_sem(r_inner_product(idx_valid), dist_circ(idx_valid), d_BinEdge);



%% Half-theory (using b coeffs from sim)
m_list = -4: 4;
r1_iHk = (0: 0.1: 10)';
[Wt_BslJ_E_cpx_coeff, Corr_FT, Corr_d_iHk] =...
    func_corr_d_theory(Wt_ctr, sigma_R, rho_C, m_list, r1_iHk);
% figure; plot(r1_iHk * d_unit, Corr_d_iHk);


figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.7, 1]); hold on;
l = zeros(1, 2);
% errorbar(d_BinCenter * d_unit, r_inner_product_1_avg, r_inner_product_1_se, 'm');
l(1) = errorbar(d_BinCenter * d_unit, r_inner_product_avg, r_inner_product_se,...
    'color', 'r', 'linewidth', 2);
xlim([0 100]); ylim([0 0.1]);
set(gca, 'XTick', 0: 10: 100, 'YTick', 0: 0.01: 0.1);
%
% plot(r1_iHk * d_unit, Corr_d_iHk * 4 * pi, 'm');
l(2) = plot(r1_iHk * d_unit, Corr_d_iHk * 2 * pi, 'linewidth', 2, 'color', 'b');
%
xlabel(['Cortical distance (\mum) (', num2str(d_unit), ' \mum per grid)']);
title('Inner product', 'fontweight', 'normal')
legend(l, {'Simulation', 'Theory'});
axis square; grid on;
ax = gca; ax.FontSize = 16;
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/Corr_d_half_theory.png']);
close;


% %% Corr(d), Sim
% dtheta_deg = 15; dphi_deg = 15;
% d_BinCenter = [1, sqrt(2), 2, sqrt(5), sqrt(8), 3.5: 0.5: 10];
% idx_fit_valid = 1: length(d_BinCenter); idx_fit_valid(5: 7) = [];
% %
% load(savename, 'K_mat', 'DoG_func', 'par_DoG_C', 'dist_circ');
% [~, ~, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se,...
%     corr_d_fitpar_ffwd_sqr, ~, ~, ~, ~, ~, ~, ~, ~] =...
% func_2D_W_stat_rho0(Wt, K_mat, DoG_func, par_DoG_C,...
%     dtheta_deg, dphi_deg, dist_circ, d_BinCenter, idx_fit_valid);
% clear K_mat DoG_func par_DoG_C dist_circ
% %
% % figure; errorbar(d_BinCenter * d_unit, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se);


% %% C(d), theory v.s. sim
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.7, 1]); hold on;
% l = zeros(1, 2);
% yyaxis left;
% l(1) = plot(r1_iHk * d_unit, Corr_d_iHk, 'linewidth', 2, 'color', 'b');
% ylim([-0.1 0.601] * 5e-5); set(gca, 'YTick', [-0.5: 0.5: 3] * 1e-5);
% yyaxis right;
% l(2) = errorbar(d_BinCenter * d_unit, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se,...
%     'linewidth', 2, 'color', 'r');
% ylim([-0.1 0.6] + 0.05); set(gca, 'YTick', 0: 0.1: 0.6);
% set(gca, 'XTick', 0: 10: 100);
% xlabel(['Cortical distance (\mum) (', num2str(d_unit), ' \mum per grid)']);
% legend(l, {'Theoretical dot product of tuning', 'Simulation tuning correlation'});
% axis square; grid on;
% ax = gca; ax.FontSize = 16;
% ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'r';
% %
% pause(1);
% print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/Corr_d_sim_theory.png']);
% %
% exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/Corr_d_sim_theory.eps'],...
%     'BackgroundColor', 'none', 'ContentType', 'vector');
% close;








%% Parsimonious model -- How \rho_K \pm \rho_C becomes C(d) in that shape
% Define coordinates
x_sqrt_fine = -(N_sqrt - 1) / 2: 0.1: (N_sqrt - 1) / 2;    % unit: grid
N_sqrt_fine = length(x_sqrt_fine);
[jj, ii] = meshgrid(x_sqrt_fine);
r_grid_fine = sqrt(jj .^ 2 + ii .^ 2); clear ii jj
% k_grid_fine = r_grid_fine / N_sqrt;
% theta_grid_fine = atan2(ii, jj);
idx_ctr_fine = (N_sqrt_fine + 1) / 2;
% r_fine = x_sqrt_fine(idx_ctr_fine: end);
r_fine = diag(r_grid_fine); r_fine = r_fine(idx_ctr_fine: end)';
rho_fine = r_fine / N_sqrt; drho_fine = rho_fine(2) - rho_fine(1);


% r1, r2 are two 'centers' of rings (s is small std)
% rho_1 = 1.5; rho_2 = 5.5;
rho_1 = (rho_K - rho_C) * N_sqrt;
rho_2 = (rho_K + rho_C) * N_sqrt;
s = 0.05;
% double-ring function
y = exp(-(r_grid_fine - rho_1).^2 / (2 * s^2)) +...
    exp(-(r_grid_fine - rho_2).^2 / (2 * s^2));
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 1]);
subplot(2, 3, 1);
imagesc(y); axis square; colorbar; clim([0 1]);
xtk = linspace(1, N_sqrt_fine, 11);
set(gca, 'XTick', xtk, 'XTickLabel', x_sqrt_fine(xtk),...
    'YTick', xtk, 'YTickLabel', x_sqrt_fine(xtk));
title({'$b(\rho) := \delta(\rho - \rho_1) + \delta(\rho - \rho_2)$',...
    '$b(|k|) = \delta(|k| - (\rho_K^0 - \rho_C^0)) + \delta(|k| - (\rho_K^0 + \rho_C^0))$'},...
    'interpreter', 'latex', 'fontweight', 'normal');
xlabel('x (1 / N)'); ylabel('x (1 / N)');
%
% sqr of cross-corr in 2D space
ycrr = xcorr2(y, y);
idx_xcorr2_ctr = (N_sqrt_fine - 1) / 2 + (1: N_sqrt_fine);
ycrr = ycrr(idx_xcorr2_ctr, idx_xcorr2_ctr);% / (N_sqrt_fine);    % should be N_sqrt_fine ^ 2
ycrr_sqr = ycrr .^ 2;
%ycrr_sqr_1D = ycrr_sqr(idx_ctr_fine, idx_ctr_fine: end);
ycrr_sqr_1D = diag(ycrr_sqr); ycrr_sqr_1D = ycrr_sqr_1D(idx_ctr_fine: end)';
%
ymax = 1500;%0.05;
subplot(2, 3, 2);
imagesc(ycrr_sqr); axis square; colorbar; clim([0 ymax]);
set(gca, 'XTick', xtk, 'XTickLabel', x_sqrt_fine(xtk),...
    'YTick', xtk, 'YTickLabel', x_sqrt_fine(xtk));
title('$(b \star b)^2(\rho)$',...
    'interpreter', 'latex', 'fontweight', 'normal');
xlabel('x (1 / N)'); ylabel('x (1 / N)');
%
subplot(2, 3, 3); hold on;
plot(rho_fine * N_sqrt, ycrr_sqr_1D);
axis([0 14 0 ymax]); set(gca, 'XTick', 0: 14);
plot(2 * rho_1 * ones(1, 2), [0 ymax], 'r--');
plot((rho_2 - rho_1) * ones(1, 2), [0 ymax], 'r--');
plot((rho_2 + rho_1) * ones(1, 2), [0 ymax], 'r--');
plot(2 * rho_2 * ones(1, 2), [0 ymax], 'r--');
title({'$(b \star b)^2(\rho)$',...
    'Location: $[2\rho_1, \rho_2 - \rho_1, \rho_2 + \rho_1, 2\rho_2]$',...
    'Location: $[2(\rho_K^0 - \rho_C^0), 2\rho_C^0, 2\rho_K^0, 2(\rho_K^0 + \rho_C^0)]$'},...
    'interpreter', 'latex', 'fontweight', 'normal');
xlabel('x (1 / N)');
%
% ratio_deltas = interp1(rho_fine * N_sqrt,...
%     ycrr_sqr_1D, [2 * rho_1, rho_2 - rho_1, rho_2 + rho_1, 2 * rho_2]);
% ratio_deltas = ratio_deltas / ratio_deltas(3);
% tmp = [rho_1, 8*rho_1*rho_2/(rho_2 - rho_1), 8*rho_1*rho_2/(rho_2 + rho_1), rho_2];
% tmp = tmp / tmp(3);


% inverse Hankel transform -- use this for theory vs sim
idx_kept = find(rho_fine <= 10 / N_sqrt);
ycrr_sqr_iHk = sum(2 * pi * drho_fine * bsxfun(@times,...
    rho_fine(idx_kept)' .* ycrr_sqr_1D(idx_kept)',...
    besselj(0, 2 * pi * rho_fine(idx_kept)' * r_fine)), 1);
%
Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;
IniVal = [1, 10, 0]; par_lb = [1e-5, 1e-5, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none');
x_fit = 0: 0.1: 100;
Err_theory = @(par) Exp2(r_fine * d_unit, par(1), par(2), par(3))...
    - ycrr_sqr_iHk;
[par_theory, ~, residual_theory, ~, ~, ~, Jacobian] =...
    lsqnonlin(Err_theory, IniVal, par_lb, par_ub, options);
CI_theory = nlparci(par_theory, residual_theory, 'jacobian', Jacobian);
par_err_theory = (CI_theory(:, 2) - CI_theory(:, 1)) / 2;
%
idx_fit_valid = 1: length(d_BinCenter); idx_fit_valid(5: 7) = [];
Err_sim = @(par) (Exp2(d_BinCenter(idx_fit_valid) * d_unit, par(1), par(2), par(3)) -...
    r_inner_product_avg(idx_fit_valid)') ./ r_inner_product_se(idx_fit_valid)';
[par_sim, ~, residual_sim, ~, ~, ~, Jacobian] =...
    lsqnonlin(Err_sim, IniVal, par_lb, par_ub, options);
CI_sim = nlparci(par_sim, residual_sim, 'jacobian', Jacobian);
par_err_sim = (CI_sim(:, 2) - CI_sim(:, 1)) / 2;
%
subplot(2, 3, 4); hold on;
l = zeros(1, 2);
yyaxis left;
l(1) = plot(r_fine * d_unit, ycrr_sqr_iHk, 'linewidth', 2, 'color', 'b');
plot(x_fit, Exp2(x_fit, par_theory(1), par_theory(2), par_theory(3)),...
    'linewidth', 1, 'color', 'b', 'linestyle', '--');
% ylim([0, 3] * 1e-3);
% ylim([0, 2.5] * 1e-3); set(gca, 'YTick', [0: 0.5: 2.5] * 1e-3);
ylim([0, 100]); set(gca, 'YTick', [0: 10: 100]);
%
ylabel('$\mathcal{F}^{-1}[(b \star b)^2](r)$', 'interpreter', 'latex');
%
yyaxis right;
l(2) = errorbar(d_BinCenter * d_unit, r_inner_product_avg, r_inner_product_se,...
    'linewidth', 2, 'color', 'r');
plot(x_fit, Exp2(x_fit, par_sim(1), par_sim(2), par_sim(3)),...
    'linewidth', 1, 'color', 'r', 'linestyle', '--');
% ylim([0 0.15]);
% ylim([0.035 0.1]); set(gca, 'YTick', [0.035, 0.04: 0.01: 0.1]);
ylim([0.03 0.1]); set(gca, 'YTick', 0.03: 0.01: 0.1);
set(gca, 'XTick', 0: 10: 100);
xlabel(['Cortical distance (\mum) (', num2str(d_unit), ' \mum per grid)']);
legend(l, {['Theory, \lambda = ', num2str(par_theory(2), '%.2f'),...
    ' \pm ', num2str(par_err_theory(2), '%.2f'), ' (\mum)'],...
    ['Simulation, \lambda = ', num2str(par_sim(2), '%.2f'),...
    ' \pm ', num2str(par_err_sim(2), '%.2f'), ' (\mum)']});
%
xlim([0 100]); set(gca, 'XTick', 0: 10: 100);
axis square; title('Inner product of tuning curves', 'fontweight', 'normal');
ax = gca; ax.FontSize = 12;
ax.YAxis(1).Color = 'b'; ax.YAxis(2).Color = 'r';

% Even most parsimonious
J0r = @(rho) besselj(0, 2 * pi * rho * r_fine / N_sqrt);
for plot_idx = [5 6]
subplot(2, 3, plot_idx); hold on;
clear l lgdtxt
l = zeros(1, 4); lgdtxt = cell(1, 4);
c_minus = J0r(rho_2 - rho_1);
l(1) = plot(r_fine * d_unit, c_minus, 'r');
c_plus = J0r(rho_2 + rho_1);
l(2) = plot(r_fine * d_unit, c_plus, 'b');
c_double = J0r(2 * rho_2);
l(3) = plot(r_fine * d_unit, (rho_2 / (4 * rho_1)) * c_double, 'm');
if plot_idx == 5
    ytmp1 = c_minus + c_plus;
    l(4) = plot(r_fine * d_unit, ytmp1 / ytmp1(1), 'k');
elseif plot_idx == 6
    ytmp2 = c_minus + c_plus + (rho_2 / (4 * rho_1)) * c_double;
    l(4) = plot(r_fine * d_unit, ytmp2 / ytmp2(1), 'k');
end
% xlim([0 14]); set(gca, 'XTick', 0: 14);
xlim([0 100]); set(gca, 'XTick', 0: 10: 100);
legend(l, {'$(\rho_K^0 + \rho_C^0) J_0(2 \pi~2 \rho_C^0~r)$',...
    '$4 (\rho_K^0 - \rho_C^0) J_0(2 \pi~2 \rho_K^0~r)$',...
    '$4 (\rho_K^0 - \rho_C^0) J_0(2 \pi~2 (\rho_K^0 + \rho_C^0) r)$',...
    'sum'}, 'location', 'northeast', 'interpreter', 'latex');
end


pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/Corr_d_parsimonious_theory.png']);
%
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/Corr_d_parsimonious_theory.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;
