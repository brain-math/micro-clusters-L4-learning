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



%% Simulation of 'canonical' parameter
filename = ['N_sqrt_', num2str(N_sqrt),...
    '_rhodeg', num2str(rho_deg), '_sA', num2str(sigma_A_input_deg), '_aC',...
    num2str(alpha_width_C), '_aK', num2str(alpha_width_K)];
savename = [dir_2D, '/results/Results_', filename, '.mat'];
N_step = 8000;
idx_trace = 0: 200: N_step;
%
% sim_2D_func(N_sqrt, rho_deg, alpha_width_C, alpha_width_K, sigma_A_input_deg,...
%     N_step, idx_trace, savename, ifUse73Save, 1);

load(savename);



%% Visualize K and C func in real and Fourier spaces
xtk_2d = 1: 5: N_sqrt; N_eig_plot = [16 48]; ylimC = [9.5 12.5]; ylimK = [0.7 0.81];
xtk_1d = 0: 2: N_sqrt - 1; ftsz_tx = 6;

figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 1]);
%
subplot(3, 4, 1);
imagesc(C_func);
cbh = colorbar; clim([-0.2 1]); set(cbh, 'YTick', -0.2: 0.2: 1);
xlabel('$\alpha_1 - \alpha_1''$', 'interpreter', 'latex');
ylabel('$\alpha_2 - \alpha_2''$', 'interpreter', 'latex');
title('$C(\vec{\alpha} - \vec{\alpha}'')$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 2);
imagesc(C_func_n);
cbh = colorbar; clim([-2 12]); set(cbh, 'YTick', -2: 2: 12);
xlabel('$n_1~(\times 1 / N)$', 'interpreter', 'latex');
ylabel('$n_2~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{C}(\vec{n})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 3); hold on;
plot(x_sqrt, C_func(idx_ctr, :));
plot([0 N_sqrt], [0 0], 'k--');
axis([0 (N_sqrt - 1) / 2 -0.2 1]); set(gca, 'YTick', -0.2: 0.2: 1);
xlabel('$|\vec{\alpha} - \vec{\alpha}''|$', 'interpreter', 'latex');
title('$C(|\vec{\alpha} - \vec{\alpha}''|)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 4); hold on;
y = C_func_n(idx_ctr, :); plot(x_sqrt, y, '-o', 'MarkerFaceColor', 'red', 'MarkerSize', 5);
plot([0 N_sqrt], [0 0], 'k--');
imax = find(y >= (max(y) - 1e-4)); plot(x_sqrt(imax(2)) * [1 1], [-2 12], 'r--');
axis([0 (N_sqrt - 1) / 2 -2 12]); set(gca, 'YTick', -2: 2: 12);
xlabel('$\rho~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{C}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
% 
%
subplot(3, 4, 5);
imagesc(K0_func);
cbh = colorbar; clim([-0.05 0.35]); set(cbh, 'YTick', [-0.05, 0, 0.05: 0.1: 0.35]);
xlabel('$x_1 - x_1''$', 'interpreter', 'latex');
ylabel('$x_2 - x_2''$', 'interpreter', 'latex');
title('$K_0(\vec{x} - \vec{x}'')$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 6);
imagesc(K0_func_n);
cbh = colorbar; clim([-0.2 0.8]); set(cbh, 'YTick', -0.2: 0.2: 0.8);
xlabel('$k_1~(\times 1 / N)$', 'interpreter', 'latex');
ylabel('$k_2~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{K_0}(\vec{k})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 7); hold on;
plot(x_sqrt, K0_func(idx_ctr, :));
plot([0 N_sqrt], [0 0], 'k--');
axis([0 (N_sqrt - 1) / 2 -0.05 0.35]); set(gca, 'YTick', [-0.05, 0, 0.05: 0.1: 0.35]);
xlabel('$|\vec{x} - \vec{x}''|$', 'interpreter', 'latex');
title('$K_0(|\vec{x} - \vec{x}''|)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 8); hold on;
y = K0_func_n(idx_ctr, :); plot(x_sqrt, y, '-o', 'MarkerFaceColor', 'red', 'MarkerSize', 5);
plot([0 N_sqrt], [0 0], 'k--');
imax = find(y >= (max(y) - 1e-4)); plot(x_sqrt(imax(2)) * [1 1], [-0.2 0.8], 'r--');
axis([0 (N_sqrt - 1) / 2 -0.2 0.8]); set(gca, 'YTick', -0.2: 0.2: 0.8);
xlabel('$d~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{K_0}(d)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
for k = 1: 8
    subplot(3, 4, k);
    axis square;
    if ismember(k, [1 2 5 6])
        set(gca, 'XTick', xtk_2d, 'XTickLabel', x_sqrt(xtk_2d),...
            'YTick', xtk_2d, 'YTickLabel', x_sqrt(xtk_2d));
    elseif ismember(k, [3 4 7 8])
        set(gca, 'XTick', xtk_1d);
    end
end
%
%
subplot(3, 4, [9 10]);
C_func_n_sorted = sortrows([j_grid(:), i_grid(:), C_func_n(:)], 3, 'descend');
plot(1: N_eig_plot(1), C_func_n_sorted(1: N_eig_plot(1), 3),...
    'color', 'k', 'marker', '.', 'markersize', 10);
for k = 1: N_eig_plot(1)
t = text(k - 0.05, C_func_n_sorted(k, 3) + 0.01,...
    ['(', num2str(C_func_n_sorted(k, 1)), ', ',...
    num2str(C_func_n_sorted(k, 2)), ')'], 'fontsize', ftsz_tx);
t.Rotation = 90;
end
axis([0 N_eig_plot(1) + 0.5 ylimC(1) ylimC(2)]); set(gca, 'XTick', 1: N_eig_plot(1));
xlabel('Index');
title('Rank~of~$\widetilde{C}(\vec{n})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, [11 12]);
K0_func_n_sorted = sortrows([j_grid(:), i_grid(:), K0_func_n(:)], 3, 'descend');
plot(1: N_eig_plot(2), K0_func_n_sorted(1: N_eig_plot(2), 3),...
    'color', 'k', 'marker', '.', 'markersize', 10);
for k = 1: N_eig_plot(2)
t = text(k - 0.05, K0_func_n_sorted(k, 3) + 0.0015,...
    ['(', num2str(K0_func_n_sorted(k, 1)), ', ',...
    num2str(K0_func_n_sorted(k, 2)), ')'], 'fontsize', ftsz_tx);
t.Rotation = 90;
end
axis([0 N_eig_plot(2) + 0.5 ylimK(1) ylimK(2)]); set(gca, 'XTick', 1: N_eig_plot(2));
xlabel('Index');
title('Rank~of~$\widetilde{K_0}(\vec{k})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/KC_func.png']);
close;
%
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/KC_func.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');


%% K and C func, finer version
dx_cont = 0.2;
x_cont = -(N_sqrt - 1) / 2: dx_cont: (N_sqrt - 1) / 2;
n_cont = x_cont / N_sqrt;
[jj, ii] = meshgrid(x_cont); r_cont = sqrt(jj .^ 2 + ii .^ 2); clear jj ii
nr_cont = r_cont / N_sqrt;
N_cont = length(x_cont); idx_ctr_cont = (N_cont + 1) / 2;
%
C_func_cont = DoG_func(r_cont, par_DoG_C);
% par(1) * exp(- d .^ 2 / (2 * par(3) ^ 2)) - par(2) * exp(- d .^ 2 / (2 * par(4) ^ 2));
DoG_func_FT = @(n, par)...
    (2 * pi * par(1) * (par(3) ^ 2)) * exp(-(2 * pi * par(3) * n) .^ 2 / 2) -...
    (2 * pi * par(2) * (par(4) ^ 2)) * exp(-(2 * pi * par(4) * n) .^ 2 / 2);
C_func_n_cont = DoG_func_FT(nr_cont, par_DoG_C);
%
W_func_n_cont = A_rec(1) * Gn_double(nr_cont, kappa(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(nr_cont, kappa(2), sigma_grid(2, 1), sigma_grid(2, 2));
K_func_n_cont = 1 ./ (1 - W_func_n_cont);
K0_func_n_cont = K_func_n_cont - 1;
%
n_cont_long = (-(N_sqrt - 1): dx_cont: (N_sqrt - 1)) / N_sqrt;
[jj, ii] = meshgrid(n_cont_long); nr_cont_long = sqrt(jj .^ 2 + ii .^ 2); clear jj ii
W_func_n_cont_long =...
    A_rec(1) * Gn_double(nr_cont_long, kappa(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(nr_cont_long, kappa(2), sigma_grid(2, 1), sigma_grid(2, 2));
K0_func_n_cont_long = (1 ./ (1 - W_func_n_cont_long)) - 1;
[r_cont_uniq, ~, r_cont_uniq_idx] = unique(r_cont);
K0_func_cont_uniq = 2 * pi * ((dx_cont / N_sqrt) ^ 2) * sum(...
    bsxfun(@times, nr_cont_long(:) .* K0_func_n_cont_long(:),...
    besselj(0, 2 * pi * nr_cont_long(:) * r_cont_uniq')), 1);
clear n_cont_long W_func_n_cont_long K0_func_n_cont_long
K0_func_cont = NaN(size(r_cont));
for i = 1: length(r_cont_uniq_idx), K0_func_cont(i) = K0_func_cont_uniq(r_cont_uniq_idx(i)); end
clear i r_cont_uniq r_cont_uniq_idx
%
A_func_cont = exp(- r_cont .^ 2 / (2 * sigma_A_input ^ 2));
A_func_n_cont = 2 * pi * (sigma_A_input ^ 2) * exp(- 2 * pi^2 * sigma_A_input^2 * nr_cont.^2);
%
n1_sup_fine = (0: 0.0001: (N_sqrt - 1) / 2) / N_sqrt;
C1_func_n_sup_fine = DoG_func_FT(n1_sup_fine, par_DoG_C);
W1_func_n_sup_fine = A_rec(1) * Gn_double(n1_sup_fine, kappa(1), sigma_grid(1, 1), sigma_grid(1, 2))...
    + A_rec(2) * Gn_double(n1_sup_fine, kappa(2), sigma_grid(2, 1), sigma_grid(2, 2));
K1_func_n_sup_fine = (1 ./ (1 - W1_func_n_sup_fine)) - 1;
A1_func_n_sup_fine = 2 * pi * (sigma_A_input ^ 2) * exp(- 2 * pi^2 * sigma_A_input^2 * n1_sup_fine.^2);


xtk_2d = linspace(1, N_cont, 5);
xtk_1d = 0: 2: N_sqrt - 1; ftsz_tx = 6;

figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 1]);
%
subplot(3, 4, 1);
imagesc(C_func_cont);
cbh = colorbar; clim([-0.2 1]); set(cbh, 'YTick', -0.2: 0.2: 1);
xlabel('$\alpha_1 - \alpha_1''$', 'interpreter', 'latex');
ylabel('$\alpha_2 - \alpha_2''$', 'interpreter', 'latex');
title('$C(\vec{\alpha} - \vec{\alpha}'')$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 2);
imagesc(C_func_n_cont);
cbh = colorbar; clim([-2 12]); set(cbh, 'YTick', -2: 2: 12);
xlabel('$n_1~(\times 1 / N)$', 'interpreter', 'latex');
ylabel('$n_2~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{C}(\vec{n})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 3); hold on;
plot(x_cont, C_func_cont(idx_ctr_cont, :));
plot([0 N_sqrt], [0 0], 'k--');
axis([0 (N_sqrt - 1) / 2 -0.2 1]); set(gca, 'YTick', -0.2: 0.2: 1);
xlabel('$|\vec{\alpha} - \vec{\alpha}''|$', 'interpreter', 'latex');
title('$C(|\vec{\alpha} - \vec{\alpha}''|)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 4); hold on;
plot(n1_sup_fine * N_sqrt, C1_func_n_sup_fine);
plot([0 N_sqrt], [0 0], 'k--');
[C1_func_n_max, idx_max] = max(C1_func_n_sup_fine); nN_max_C = n1_sup_fine(idx_max) * N_sqrt;
plot(nN_max_C * [1 1], [-2 12], 'r--');
axis([0 (N_sqrt - 1) / 2 -2 12]); set(gca, 'YTick', -2: 2: 12);
xlabel({'$\rho~(\times 1 / N)$', ['$\rho_C^0 = ', num2str(nN_max_C, '%.3f'),...
    '/N, \widetilde{C}_{\rm max.} = ', num2str(C1_func_n_max, '%.3f'), '$']}, 'interpreter', 'latex');
title('$\widetilde{C}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
% 
%
subplot(3, 4, 5);
imagesc(K0_func_cont);
cbh = colorbar; clim([-0.1 1.05]); set(cbh, 'YTick', [-0.1, 0: 0.2: 1]);
xlabel('$x_1 - x_1''$', 'interpreter', 'latex');
ylabel('$x_2 - x_2''$', 'interpreter', 'latex');
title('$K_0(\vec{x} - \vec{x}'')$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 6);
imagesc(K0_func_n_cont);
cbh = colorbar; clim([-0.2 0.8]); set(cbh, 'YTick', -0.2: 0.2: 0.8);
xlabel('$k_1~(\times 1 / N)$', 'interpreter', 'latex');
ylabel('$k_2~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{K_0}(\vec{k})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 7); hold on;
plot(x_cont, K0_func_cont(idx_ctr_cont, :));
plot([0 N_sqrt], [0 0], 'k--');
axis([0 (N_sqrt - 1) / 2 -0.1 1.05]); set(gca, 'YTick', [-0.1, 0: 0.2: 1]);
xlabel('$|\vec{x} - \vec{x}''|$', 'interpreter', 'latex');
title('$K_0(|\vec{x} - \vec{x}''|)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 8); hold on;
plot(n1_sup_fine * N_sqrt, K1_func_n_sup_fine);
plot([0 N_sqrt], [0 0], 'k--');
[K1_func_n_max, idx_max] = max(K1_func_n_sup_fine); nN_max_K = n1_sup_fine(idx_max) * N_sqrt;
plot(nN_max_K * [1 1], [-0.6 0.8], 'r--');
axis([0 (N_sqrt - 1) / 2 -0.6 0.8]); set(gca, 'YTick', -0.6: 0.2: 0.8);
xlabel({'$d~(\times 1 / N)$', ['$d_K^0 = ', num2str(nN_max_K, '%.3f'),...
    '/N, \widetilde{K_0}_{\rm max.} = ',...
    num2str(K1_func_n_max, '%.3f'), '$']}, 'interpreter', 'latex');
title('$\widetilde{K_0}(d)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 9);
imagesc(A_func_cont);
cbh = colorbar; clim([0 1.01]); set(cbh, 'YTick', 0: 0.2: 1);
xlabel('$x_1 - \alpha_1''$', 'interpreter', 'latex');
ylabel('$x_2 - \alpha_2''$', 'interpreter', 'latex');
title('$A(\vec{x} - \vec{\alpha})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 10);
imagesc(A_func_n_cont);
cbh = colorbar; clim([0 300]); set(cbh, 'YTick', 0: 50: 300);
xlabel('$k_1~(\times 1 / N)$', 'interpreter', 'latex');
ylabel('$k_2~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{A}(\vec{n})$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 11); hold on;
A_func_cont_long = diag(A_func_cont)';
plot(x_cont * sqrt(2), A_func_cont_long);
axis([0 14 0 1.01]); set(gca, 'XTick', 0: 2: 14, 'YTick', 0: 0.2: 1);
xlabel('$|\vec{x} - \vec{\alpha}|$', 'interpreter', 'latex');
title('$A(|\vec{x} - \vec{\alpha}|)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 4, 12); hold on;
plot(n1_sup_fine * N_sqrt, A1_func_n_sup_fine);
axis([0 (N_sqrt - 1) / 2 0 310]); set(gca, 'YTick', 0: 50: 300);
xlabel('$\rho~(\times 1 / N)$', 'interpreter', 'latex');
title('$\widetilde{A}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
for k = 1: 12
    subplot(3, 4, k);
    axis square;
    if ismember(k, [1 2 5 6 9 10])
        set(gca, 'XTick', xtk_2d, 'XTickLabel', x_cont(xtk_2d),...
            'YTick', xtk_2d, 'YTickLabel', x_cont(xtk_2d));
    elseif ismember(k, [3 4 7 8 12])
        set(gca, 'XTick', xtk_1d);
    end
end
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/KCA_func_fine.png']);
%
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/KCA_func_fine.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;



d_unit = (rho_mag * rho_deg);
xmax = 50;    % 1000 / 21
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 1]);
%
subplot(3, 3, 1); hold on;
plot(n1_sup_fine * (1000 / d_unit), C1_func_n_sup_fine);
plot([0 xmax], [0 0], 'k--');
[C1_func_n_max, idx_max] = max(C1_func_n_sup_fine); nN_max_C = n1_sup_fine(idx_max) * (1000 / d_unit);
plot(nN_max_C * [1 1], [-2 12], 'r--');
axis square;
axis([0 xmax -2 12]); set(gca, 'XTick', 0: 5: xmax, 'YTick', -2: 2: 12);
xlabel({'$\rho (mm^{-1})$', ['$\rho_C^0 = ', num2str(nN_max_C, '%.3f'),...
    ' (mm^{-1}), \widetilde{C}_{\rm max.} = ',...
    num2str(C1_func_n_max, '%.3f'), '$']}, 'interpreter', 'latex');
title('$\widetilde{C}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 3, 4); hold on;
plot(n1_sup_fine * (1000 / d_unit), K1_func_n_sup_fine + 1);
plot([0 xmax], [1 1], 'k--');
[K1_func_n_max, idx_max] = max(K1_func_n_sup_fine); nN_max_K = n1_sup_fine(idx_max) * (1000 / d_unit);
plot(nN_max_K * [1 1], [0.4 1.8], 'r--');
axis square;
axis([0 xmax 0.4 1.8]); set(gca, 'XTick', 0: 5: xmax, 'YTick', 0.4: 0.2: 1.8);
xlabel({'$\rho (mm^{-1})$', ['$\rho_K^0 = ', num2str(nN_max_K, '%.3f'),...
    ' (mm^{-1}) \widetilde{K}_{\rm max.} = ',...
    num2str(K1_func_n_max + 1, '%.3f'), '$']}, 'interpreter', 'latex');
title('$\widetilde{K}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 3, 7); hold on;
plot(n1_sup_fine * (1000 / d_unit), A1_func_n_sup_fine);
axis square;
axis([0 xmax 0 310]); set(gca, 'XTick', 0: 5: xmax, 'YTick', 0: 50: 300);
xlabel('$\rho (mm^{-1})$', 'interpreter', 'latex');
title('$\widetilde{A}(\rho)$', 'FontWeight', 'normal', 'interpreter', 'latex');
%
subplot(3, 3, [2 3 5 6 8 9]); hold on;
plot(n1_sup_fine * (1000 / d_unit), C1_func_n_sup_fine, 'color', [0 0.5 1]);
plot(n1_sup_fine * (1000 / d_unit), K1_func_n_sup_fine + 1, 'color', [1 0.5 0]);
plot(nN_max_C * [1 1], [-2 12], 'color', [0 0.5 1], 'linestyle', '--');
plot(nN_max_K * [1 1], [-2 12], 'color', [1 0.5 0], 'linestyle', '--');
axis square;
axis([0 xmax -2 12]); set(gca, 'XTick', 0: 5: xmax, 'YTick', -2: 2: 12);
xlabel('$\rho (mm^{-1})$', 'interpreter', 'latex');
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/KCA_FT.png']);
%
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/KCA_FT.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;



%% Visualize RF over t -- also absFT_sqr
d_boundary = 2; cmax = 0.015; cmax_ftsqr = 0.25;

dir_RF = [dir_2D, '/figures/', filename, '/RF/RF'];
dir_RF_absFTsqr = [dir_2D, '/figures/', filename, '/RF/RF_absFTsqr'];
if ~exist(dir_RF, 'dir'), mkdir(dir_RF); end
if ~exist(dir_RF_absFTsqr, 'dir'), mkdir(dir_RF_absFTsqr); end

W_trace_ctr = NaN(size(W_trace));
for ti = 1: N_trace
    W_trace_ctr(:, :, ti) = circmat_align_2D(W_trace(:, :, ti), 1);
end
clear ti
%
N_sqrt_cut_half = 5; N_sqrt_cut = 2 * N_sqrt_cut_half + 1;
idx_cut = (N_sqrt - N_sqrt_cut) / 2 + [1: N_sqrt_cut];
W_ctr_absFTsqr_trace = NaN(N, N_sqrt_cut ^ 2, N_trace);
for ti = 1: N_trace
for xi = 1: N
    Wx = reshape(W_trace_ctr(xi, :, ti), [N_sqrt, N_sqrt]);
    ft_tmp = abs(fftshift(fft2(ifftshift(Wx))));
    ft_tmp = ft_tmp(idx_cut, idx_cut);
    W_ctr_absFTsqr_trace(xi, :, ti) = (ft_tmp(:) .^ 2)';
end
end
clear ti xi Wx ft_tmp idx_cut

for plot_i = 1: N_trace
    RF_2D_visualize_fast(W_trace_ctr(:, :, plot_i), d_boundary, -cmax, cmax, Clr_RWB);
    %
    pause(0.5);
    print(gcf, '-dpng', [dir_RF, '/RF_', num2str(idx_trace(plot_i)), '.png']);
    close;
end
%
for plot_i = 1: N_trace
    RF_2D_visualize_fast(W_ctr_absFTsqr_trace(:, :, plot_i), d_boundary, 0, cmax_ftsqr, 0);
    %
    pause(0.5);
    print(gcf, '-dpng', [dir_RF_absFTsqr, '/RF_absFTsqr_', num2str(idx_trace(plot_i)), '.png']);
    close;
end
%
% video for RF(t)
v = VideoWriter(fullfile([dir_2D, '/figures/', filename, '/RF'], 'RF.mp4'), 'MPEG-4');
v.FrameRate = 8;
open(v);
for plot_i = 1: N_trace
    im0 = imread([dir_RF, '/RF_', num2str(idx_trace(plot_i)), '.png']);
    im = insertText(im0, [10, 0], ['Frame ', num2str(idx_trace(plot_i))], ...
        'FontSize', 20, 'BoxOpacity', 0.4, ...
        'BoxColor', 'black', 'TextColor', 'white');
    % figure; imshow(im)
    writeVideo(v, im2frame(im));
end
writeVideo(v, im2frame(im));
close(v);



%% Output statistics
trace_idx_select = find(ismember(idx_trace, [1000, 2000, 4000, 8000]));
%
dtheta_deg = 15; dphi_deg = 15;
d_BinCenter = [1, sqrt(2), 2, sqrt(5), sqrt(8), 3.5: 0.5: 10];
%
idx_fit_valid = 1: length(d_BinCenter); idx_fit_valid(5: 7) = [];
%
% Still use optimal rho_0 as the only 
[nC_max_fine, r_ffwd_sqr, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se,...
    corr_d_fitpar_ffwd_sqr, theta_pref_ffwd_sqr, gOSI_ffwd_sqr,...
    r_rec, corr_d_rec_avg, corr_d_rec_se,...
    corr_d_fitpar_rec, theta_pref_rec, gOSI_rec] =...
    func_2D_W_stat_rho0(W_trace(:, :, trace_idx_select),...
    K_mat, DoG_func, par_DoG_C, dtheta_deg, dphi_deg, dist_circ, d_BinCenter, idx_fit_valid);
%
% n0_list = [1: 0.5: 4] / N_sqrt;
% thr_valid_corr = 0.05;
% [r_ffwd_sqr, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se,...
%     corr_d_fitpar_ffwd_sqr, theta_pref_ffwd_sqr, gOSI_ffwd_sqr,...
%     r_rec, corr_d_rec_avg, corr_d_rec_se,...
%     corr_d_fitpar_rec, theta_pref_rec, gOSI_rec] =...
%     func_2D_W_stat_n_all(W_trace(:, :, trace_idx_select),...
%     K_mat, n0_list, dtheta_deg, dphi_deg, dist_circ, d_BinCenter, idx_fit_valid, thr_valid_corr);


t_i_mid = find(trace_idx_select == find(idx_trace == 1000));
t_i_end = length(trace_idx_select);
%
theta_stim_deg_ext = 0: dtheta_deg: 180;
N_curve_plot = 50; idx_plot = randperm(N); idx_plot = sort(idx_plot(1: N_curve_plot));
%
d_unit = (rho_mag * rho_deg);
N_idx_select = length(trace_idx_select); clr = jet(N_idx_select) * 0.875;
d_fit = 0: 1: 100; Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;
%
for plot_i = 1: 2
%
if plot_i == 1
    rsp = r_ffwd_sqr;
    theta_pref = theta_pref_ffwd_sqr;
    gOSI = gOSI_ffwd_sqr;
    corr_d_avg = corr_d_ffwd_sqr_avg;
    corr_d_se = corr_d_ffwd_sqr_se;
    corr_d_fitpar = corr_d_fitpar_ffwd_sqr;
    filename_append = 'ffwd_sqr';
elseif plot_i == 2
    rsp = r_rec;
    theta_pref = theta_pref_rec;
    gOSI = gOSI_rec;
    corr_d_avg = corr_d_rec_avg;
    corr_d_se = corr_d_rec_se;
    corr_d_fitpar = corr_d_fitpar_rec;
    filename_append = 'rec';
end
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.8, 1]);
%
% Tuning curve bunches
subplot(3, 3, 1); hold on;
for k = 1: N_curve_plot
    y = rsp(idx_plot(k), :, t_i_mid);
    plot(theta_stim_deg_ext, [y(end) y]);
end
axis([0, 180, 0, ceil(max(max(rsp(idx_plot, :, t_i_mid))) / 0.05) * 0.05]);
set(gca, 'XTick', 0: 15: 180);
xlabel('Stimulus orientation (deg.)');
ylabel(['Step ', num2str(idx_trace(trace_idx_select(t_i_mid)))]);
title('Tuning curves', 'fontweight', 'normal');
%
subplot(3, 3, 4); hold on;
for k = 1: N_curve_plot
    y = rsp(idx_plot(k), :, t_i_end);
    plot(theta_stim_deg_ext, [y(end) y]);
end
axis([0, 180, 0, ceil(max(max(rsp(idx_plot, :, t_i_end))) / 0.05) * 0.05]);
set(gca, 'XTick', 0: 15: 180);
xlabel('Stimulus orientation (deg.)');
ylabel(['Step ', num2str(idx_trace(trace_idx_select(t_i_end)))]);
title('Tuning curves', 'fontweight', 'normal');
%
%
% Orientation map * gOSI as alpha value
subplot(3, 3, 2); hold on; box on;
h = imagesc(reshape(theta_pref(:, t_i_mid), [N_sqrt, N_sqrt]));
set(h, 'alphadata', reshape(gOSI(:, t_i_mid), [N_sqrt, N_sqrt])); 
title('Orientation map (gOSI as transparency)', 'fontweight', 'normal');
%
subplot(3, 3, 5); hold on; box on;
h = imagesc(reshape(theta_pref(:, t_i_end), [N_sqrt, N_sqrt]));
set(h, 'alphadata', reshape(gOSI(:, t_i_end), [N_sqrt, N_sqrt]));
title('Orientation map (gOSI as transparency)', 'fontweight', 'normal');
%
for subplot_k = [2 5]
    subplot(3, 3, subplot_k);
    axis square; axis([0.5, N_sqrt + 0.5, 0.5, N_sqrt + 0.5]);
    colormap(gca, 'hsv');
    clb = colorbar; clim([-1e-4, pi + 1e-4]); clb.Ticks = 0: pi/4: pi;
    clb.TickLabels = {'0', '45^o', '90^o', '135^o', '180^o'};
    set(gca, 'XTick', [1, N_sqrt], 'YTick', [1, N_sqrt]);
end
%
subplot(6, 3, 3); hold on;
histogram(theta_pref(:, t_i_mid) * (180 / pi), 0: dtheta_deg: 180);
xlim([0, 180]); set(gca, 'XTick', 0: dtheta_deg: 180); ylabel({'Pref.', 'Orientation'});
subplot(6, 3, 6); hold on;
histogram(gOSI(:, t_i_mid), 0: 0.05: 1);
xlim([0, 1]); set(gca, 'XTick', 0: 0.1: 1); ylabel('gOSI');
subplot(6, 3, 9); hold on;
histogram(theta_pref(:, t_i_end) * (180 / pi), 0: dtheta_deg: 180);
xlim([0, 180]); set(gca, 'XTick', 0: dtheta_deg: 180); ylabel({'Pref.', 'Orientation'});
subplot(6, 3, 12); hold on;
histogram(gOSI(:, t_i_end), 0: 0.05: 1);
xlim([0, 1]); set(gca, 'XTick', 0: 0.1: 1); ylabel('gOSI');
%
% 
% Corr(d)
subplot(3, 3, 9); hold on;
l = zeros(1, N_idx_select); lgdtxt = cell(1, N_idx_select);
for k = 1: N_idx_select
    l(k) = errorbar(d_BinCenter * d_unit,...
        corr_d_avg(:, k), corr_d_se(:, k), 'color', clr(k, :), 'linewidth', 1);
    if (corr_d_avg(1, k) >= 0.15)
        par = corr_d_fitpar(:, k, 1); par(2) = par(2) * d_unit;
        par_err = corr_d_fitpar(:, k, 2); par_err(2) = par_err(2) * d_unit;
        plot(d_fit, Exp2(d_fit, par(1), par(2), par(3)),...
            'linestyle', '--', 'color', clr(k, :), 'linewidth', 0.5);
        lgdtxt{k} = ['Step ', num2str(idx_trace(trace_idx_select(k))), ', A = ',...
            num2str(par(1), '%.3f'), ' \pm ', num2str(par_err(1), '%.3f'), ', \lambda = ',...
            num2str(par(2), '%.2f'), ' \pm ', num2str(par_err(2), '%.3f'), ' (\mum)'];
    else
        lgdtxt{k} = ['Step ', num2str(idx_trace(trace_idx_select(k)))];
    end
end
grid on;
plot([0 100], [0 0], 'k--');
axis([0 100 -0.1 0.8]); set(gca, 'XTick', 0: 10: 100, 'YTick', [-0.1 0: 0.2: 0.8]);
xlabel(['Cortical distance (\mum) (', num2str(d_unit), ' \mum per grid)']);
title({'Tuning curve correlation', ['(max. rsp. n = ',...
    num2str(nC_max_fine * N_sqrt, '%.2f'),...
    ' / ', num2str(N_sqrt), ' grids)']}, 'fontweight', 'normal');
hlgd = legend(l, lgdtxt);
set(hlgd, 'Position', [0.35, 0.1, 0.3, 0.08]);    % [left bottom width height]
%
pause(2);
print(gcf, '-dpng', [dir_2D, '/figures/', filename, '/W_stat_', filename_append, '.png']);
close;

end
%
% exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/W_stat.eps'],...
%     'BackgroundColor', 'none', 'ContentType', 'vector');
% exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/W_stat_corr.eps'],...
%     'BackgroundColor', 'none', 'ContentType', 'vector');
