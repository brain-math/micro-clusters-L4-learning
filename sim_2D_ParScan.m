%% SIM
dir_0 = 'D:/Study/CompNeuro/Projects/L4_Oja';
dir_2D = [dir_0, '/sim_theory_2D_fin'];
addpath(dir_0); addpath(dir_2D);
addpath(genpath('D:/Study/CompNeuro/Projects/Functions_simul'));

N_sqrt = 21; N = N_sqrt ^ 2;
rho_deg = 0.5;
sigma_A_input_deg_0 = 3.5;
ifUse73Save = 0;
%
N_step = 4000;
idx_trace = [0, N_step];
%
filename0 = ['N_sqrt_', num2str(N_sqrt), '_rhodeg', num2str(rho_deg),...
    '_sA', num2str(sigma_A_input_deg_0), '_aC1_aK1'];
load([dir_2D, '/results/Results_', filename0, '.mat'], 'dist_circ');
dtheta_deg = 15; N_theta = round(180 / dtheta_deg); dphi_deg = 15;
d_BinCenter = [1, sqrt(2), 2, sqrt(5), sqrt(8), 3.5: 0.5: 10];
idx_fit_valid = [1: 3, length(d_BinCenter) + (-3: 0)];
% idx_fit_valid = 1: length(d_BinCenter); idx_fit_valid(5: 7) = [];
%
N_loop = 5;
%
dir_save = [dir_2D, '/results/ParScan_N_sqrt_', num2str(N_sqrt)];



%% 1. lA fixed, lK v.s. lC
% 11 x 11 x 5, 166.13 min
% X For (aC vs aK) consider broad aK to 5?
sigma_A_input_deg = sigma_A_input_deg_0;
lA = sigma_A_input_deg / sigma_A_input_deg_0;
lC_list = 0.8: 0.1: 1.8;
lK_list = 0.8: 0.1: 1.8;
%
tic;
for i = 1: length(lC_list)
for j = 1: length(lK_list)
    lC = lC_list(i);
    lK = lK_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    r_ffwd_sqr_rounds = NaN(N, N_theta, N_loop);
    for loop_k = 1: N_loop
        savename = [dir_save, '/Results_', filename, '_round', num2str(loop_k), '.mat'];
        sim_2D_func(N_sqrt, rho_deg, lC, lK, sigma_A_input_deg,...
            N_step, idx_trace, savename, ifUse73Save, 0);
        %
        load(savename, 'DoG_func', 'par_DoG_C', 'W_trace');
        Wt = W_trace(:, :, end);
        [~, r_ffwd_sqr_rounds(:, :, loop_k)] =...
            func_2D_W_stat_simple_rsp(Wt, DoG_func, par_DoG_C, dtheta_deg, dphi_deg);
        %
        fprintf(['ParScan N = ', num2str(N_sqrt), ', Group1, (',...
            num2str(i), ', ', num2str(j), ')  / (', num2str(length(lC_list)),...
            ', ', num2str(length(lK_list)), '), loop ', num2str(loop_k), ' / ',...
            num2str(N_loop), ' fin, ', num2str(toc / 60, '%.2f'), ' min.\n']);
        delete(savename);
        clear savename DoG_func par_DoG_C W_trace Wt
    end
    %
    [corr_d_avg, corr_d_se, corr_d_fitpar, ip_d_avg, ip_d_se, ip_d_fitpar] =...
        func_2D_W_stat_simple_CorrMultip(...
        r_ffwd_sqr_rounds, dist_circ, d_BinCenter, idx_fit_valid);
    save([dir_save, '/Results_integration_', filename, '.mat'],...
        'r_ffwd_sqr_rounds', 'corr_d_avg', 'corr_d_se',...
        'corr_d_fitpar', 'ip_d_avg', 'ip_d_se', 'ip_d_fitpar');
    clear lC lK filename
end
end
clear sigma_A_input_deg lA lC_list lK_list i j


%% 2. lC fixed, lK v.s. lA
lC = 1;
lA_list = 0.5: 0.25: 2;
lK_list = 0.8: 0.1: 1.8;
%
tic;
for i = 1: length(lA_list)
for j = 1: length(lK_list)
    lA = lA_list(i);
    sigma_A_input_deg = sigma_A_input_deg_0 * lA;
    lK = lK_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    r_ffwd_sqr_rounds = NaN(N, N_theta, N_loop);
    for loop_k = 1: N_loop
        savename = [dir_save, '/Results_', filename, '_round', num2str(loop_k), '.mat'];
        sim_2D_func(N_sqrt, rho_deg, lC, lK, sigma_A_input_deg,...
            N_step, idx_trace, savename, ifUse73Save, 0);
        %
        load(savename, 'DoG_func', 'par_DoG_C', 'W_trace');
        Wt = W_trace(:, :, end);
        [~, r_ffwd_sqr_rounds(:, :, loop_k)] =...
            func_2D_W_stat_simple_rsp(Wt, DoG_func, par_DoG_C, dtheta_deg, dphi_deg);
        %
        fprintf(['ParScan N = ', num2str(N_sqrt), ', Group2, (',...
            num2str(i), ', ', num2str(j), ')  / (', num2str(length(lA_list)),...
            ', ', num2str(length(lK_list)), '), loop ', num2str(loop_k), ' / ',...
            num2str(N_loop), ' fin, ', num2str(toc / 60, '%.2f'), ' min.\n']);
        delete(savename);
        clear savename DoG_func par_DoG_C W_trace Wt
    end
    %
    [corr_d_avg, corr_d_se, corr_d_fitpar, ip_d_avg, ip_d_se, ip_d_fitpar] =...
        func_2D_W_stat_simple_CorrMultip(...
        r_ffwd_sqr_rounds, dist_circ, d_BinCenter, idx_fit_valid);
    save([dir_save, '/Results_integration_', filename, '.mat'],...
        'r_ffwd_sqr_rounds', 'corr_d_avg', 'corr_d_se',...
        'corr_d_fitpar', 'ip_d_avg', 'ip_d_se', 'ip_d_fitpar');
    clear lA sigma_A_input_deg lK filename
end
end
clear lC lA_list lK_list i j


%% 3. lK fixed, lC v.s. lA
lK = 1;
lA_list = 0.5: 0.25: 2;
lC_list = 0.8: 0.1: 1.8;
%
tic;
for i = 1: length(lA_list)
for j = 1: length(lC_list)
    lA = lA_list(i);
    sigma_A_input_deg = sigma_A_input_deg_0 * lA;
    lC = lC_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    r_ffwd_sqr_rounds = NaN(N, N_theta, N_loop);
    for loop_k = 1: N_loop
        savename = [dir_save, '/Results_', filename, '_round', num2str(loop_k), '.mat'];
        sim_2D_func(N_sqrt, rho_deg, lC, lK, sigma_A_input_deg,...
            N_step, idx_trace, savename, ifUse73Save, 0);
        %
        load(savename, 'DoG_func', 'par_DoG_C', 'W_trace');
        Wt = W_trace(:, :, end);
        [~, r_ffwd_sqr_rounds(:, :, loop_k)] =...
            func_2D_W_stat_simple_rsp(Wt, DoG_func, par_DoG_C, dtheta_deg, dphi_deg);
        %
        fprintf(['ParScan N = ', num2str(N_sqrt), ', Group2, (',...
            num2str(i), ', ', num2str(j), ')  / (', num2str(length(lA_list)),...
            ', ', num2str(length(lC_list)), '), loop ', num2str(loop_k), ' / ',...
            num2str(N_loop), ' fin, ', num2str(toc / 60, '%.2f'), ' min.\n']);
        delete(savename);
        clear savename DoG_func par_DoG_C W_trace Wt
    end
    %
    [corr_d_avg, corr_d_se, corr_d_fitpar, ip_d_avg, ip_d_se, ip_d_fitpar] =...
        func_2D_W_stat_simple_CorrMultip(...
        r_ffwd_sqr_rounds, dist_circ, d_BinCenter, idx_fit_valid);
    save([dir_save, '/Results_integration_', filename, '.mat'],...
        'r_ffwd_sqr_rounds', 'corr_d_avg', 'corr_d_se',...
        'corr_d_fitpar', 'ip_d_avg', 'ip_d_se', 'ip_d_fitpar');
    clear lC lA sigma_A_input_deg filename
end
end
clear lK lC_list lA_list i j




%% Sim plots
load([dir_2D, '/results/Results_N_sqrt_', num2str(N_sqrt), '_rhodeg',...
    num2str(rho_deg), '_sA', num2str(sigma_A_input_deg_0), '_aC1_aK1.mat'],...
    'rho_mag');
d_unit = rho_mag * rho_deg;
d_fit = 0: 1: 100; Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;


%% 1)
sigma_A_input_deg = sigma_A_input_deg_0;
lA = sigma_A_input_deg / sigma_A_input_deg_0;
lC_list = 1.8: -0.1: 0.8; N_lC = length(lC_list);
lK_list = 0.8: 0.1: 1.8; N_lK = length(lK_list);
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.8, 1]);
suptitle(['Tuning correlation (lA = 1, Max rsp. n, t = ', num2str(N_step), ')'], 2, 0.97);
%
sigma_corr = NaN(N_lC, N_lK);
for i = 1: N_lC
for j = 1: N_lK
    lC = lC_list(i);
    lK = lK_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    load([dir_save, '/Results_integration_', filename, '.mat'],...
        'corr_d_avg', 'corr_d_se', 'corr_d_fitpar');
    subplot(N_lC, N_lK, (i - 1) * N_lK + j); hold on;
    plot(d_BinCenter' * d_unit, corr_d_avg(:, end),...
        'color', 'b', 'linewidth', 1, 'marker', '.', 'markersize', 5);
    par = corr_d_fitpar(:, end, 1); par(2) = par(2) * d_unit; sigma_corr(i, j) = par(2);
    plot(d_fit, Exp2(d_fit, par(1), par(2), par(3)),...
        'linestyle', '--', 'color', [0.5 0 1], 'linewidth', 0.5);
    plot([0 100], [0 0], 'k--');
    axis([0 100 -0.1 1]); set(gca, 'XTick', 0: 10: 100, 'YTick', [-0.1 0: 0.25: 1],...
        'XTickLabel', {'0', '', '20', '', '40', '', '60', '', '80', '', '100'},...
        'YTickLabel', {'', '0', '0.25', '0.5', '0.75', '1'}); grid on;
    if j ~= 1, set(gca, 'YTickLabel', []); end
    if i ~= length(lK_list), set(gca, 'XTickLabel', []); end
    ax = gca; ax.FontSize = 6;    % for ticks
    %
    text(95, 0.85, ['\sigma = ', num2str(par(2), '%.2f')],...
        'HorizontalAlignment', 'right', 'FontSize', 7);
    text(95, 0.55, '(\mum)', 'HorizontalAlignment', 'right', 'FontSize', 7);
    if (i == N_lK) & (j == 1)
        xlabel({'Cortical distance (\mum)',...
            ['(', num2str(d_unit), ' \mum per grid)']}, 'FontSize', 8);
    end
    if j == 1, ylabel(['lC = ', num2str(lC)], 'FontSize', 10); end
    if i == 1, title(['lK = ', num2str(lK)], 'fontweight', 'normal', 'FontSize', 10); end
end
end
%
sigma_corr_sim = sigma_corr;
save([dir_2D, '/results/ParScan_N_sqrt_', num2str(N_sqrt), '/Theory_lC_lK.mat'],...
    'sigma_corr_sim', 'lC_list', 'lK_list', '-append');
%
pause(2);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lK_vs_lC_supp.png']);
close;


lK_list_fine = 0.8: 0.05: 1.8;
lC_list_fine = 1.8: -0.05: 0.8;
[x, y] = meshgrid(lK_list, lC_list);
[x_fine, y_fine] = meshgrid(lK_list_fine, lC_list_fine);
sigma_corr_fine = interp2(x, y, sigma_corr, x_fine, y_fine);
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 0.8]);
subplot(1, 2, 1); imagesc(sigma_corr);
set(gca, 'XTick', 1: N_lK, 'XTickLabel', lK_list,...
    'YTick', 1: N_lC, 'YTickLabel', lC_list);
subplot(1, 2, 2); imagesc(sigma_corr_fine);
idx = find(ismember(lK_list_fine, lK_list));
set(gca, 'XTick', idx, 'XTickLabel', lK_list_fine(idx),...
    'YTick', idx, 'YTickLabel', lC_list_fine(idx));
for k = 1: 2
    subplot(1, 2, k);
    axis square; h = colorbar;
    clim([5 13]); % clim([floor(min(sigma_corr(:))) ceil(max(sigma_corr(:)))]);
    xlabel('lK'); ylabel('lC'); title(h, '\sigma (\mum)');
end
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lK_vs_lC.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/ParScan_lK_vs_lC.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;


%% 2) -- (0.75, 1.75) rather than (0.5, 2)
lC = 1;
lA_list = 1.75: -0.25: 0.75; N_lA = length(lA_list);
lK_list = 0.8: 0.1: 1.8; N_lK = length(lK_list);
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.8, 1]);
suptitle(['Tuning correlation (lC = 1, Max rsp. n, t = ', num2str(N_step), ')'], 2, 0.97);
%
sigma_corr = NaN(N_lA, N_lK);
for i = 1: N_lA
for j = 1: N_lK
    lA = lA_list(i);
    lK = lK_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    load([dir_save, '/Results_integration_', filename, '.mat'],...
        'corr_d_avg', 'corr_d_se', 'corr_d_fitpar');
    subplot(N_lA, N_lK, (i - 1) * N_lK + j); hold on;
    plot(d_BinCenter' * d_unit, corr_d_avg(:, end),...
        'color', 'b', 'linewidth', 1, 'marker', '.', 'markersize', 5);
    par = corr_d_fitpar(:, end, 1); par(2) = par(2) * d_unit; sigma_corr(i, j) = par(2);
    plot(d_fit, Exp2(d_fit, par(1), par(2), par(3)),...
        'linestyle', '--', 'color', [0.5 0 1], 'linewidth', 0.5);
    plot([0 100], [0 0], 'k--');
    axis([0 100 -0.1 1]); set(gca, 'XTick', 0: 10: 100, 'YTick', [-0.1 0: 0.25: 1],...
        'XTickLabel', {'0', '', '20', '', '40', '', '60', '', '80', '', '100'},...
        'YTickLabel', {'', '0', '0.25', '0.5', '0.75', '1'}); grid on;
    if j ~= 1, set(gca, 'YTickLabel', []); end
    if i ~= length(lK_list), set(gca, 'XTickLabel', []); end
    ax = gca; ax.FontSize = 6;    % for ticks
    %
    text(95, 0.85, ['\sigma = ', num2str(par(2), '%.2f')],...
        'HorizontalAlignment', 'right', 'FontSize', 7);
    text(95, 0.55, '(\mum)', 'HorizontalAlignment', 'right', 'FontSize', 7);
    if (i == N_lK) & (j == 1)
        xlabel({'Cortical distance (\mum)',...
            ['(', num2str(d_unit), ' \mum per grid)']}, 'FontSize', 8);
    end
    if j == 1, ylabel(['lA = ', num2str(lA)], 'FontSize', 10); end
    if i == 1, title(['lK = ', num2str(lK)], 'fontweight', 'normal', 'FontSize', 10); end
end
end
%
pause(2);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lK_vs_lA_supp.png']);
close;


lK_list_fine = 0.8: 0.05: 1.8;
lA_list_fine = 1.75: -0.125: 0.75;
[x, y] = meshgrid(lK_list, lA_list);
[x_fine, y_fine] = meshgrid(lK_list_fine, lA_list_fine);
sigma_corr_fine = interp2(x, y, sigma_corr, x_fine, y_fine);
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 0.8]);
subplot(1, 2, 1); imagesc(sigma_corr);
set(gca, 'XTick', 1: N_lK, 'XTickLabel', lK_list,...
    'YTick', 1: N_lA, 'YTickLabel', lA_list);
subplot(1, 2, 2); imagesc(sigma_corr_fine);
idxK = find(ismember(lK_list_fine, lK_list));
idxA = find(ismember(lA_list_fine, lA_list));
set(gca, 'XTick', idxK, 'XTickLabel', lK_list_fine(idxK),...
    'YTick', idxA, 'YTickLabel', lA_list_fine(idxA));
for k = 1: 2
    subplot(1, 2, k);
    axis square; h = colorbar;
    clim([5 11]); % clim([floor(min(sigma_corr(:))) ceil(max(sigma_corr(:)))]);
    xlabel('lK'); ylabel('lA'); title(h, '\sigma (\mum)');
end
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lK_vs_lA.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/ParScan_lK_vs_lA.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;


%% 3) -- (0.75, 1.75) rather than (0.5, 2)
lK = 1;
lA_list = 1.75: -0.25: 0.75; N_lA = length(lA_list);
lC_list = 0.8: 0.1: 1.8; N_lC = length(lC_list);
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.8, 1]);
suptitle(['Tuning correlation (lK = 1, Max rsp. n, t = ', num2str(N_step), ')'], 2, 0.97);
%
sigma_corr = NaN(N_lA, N_lC);
for i = 1: N_lA
for j = 1: N_lC
    lA = lA_list(i);
    lC = lC_list(j);
    %
    filename = ['N_sqrt_', num2str(N_sqrt),...
        '_rhodeg', num2str(rho_deg, '%.1f'), '_aA', num2str(lA, '%.2f'), '_aC',...
        num2str(lC, '%.2f'), '_aK', num2str(lK, '%.2f')];
    load([dir_save, '/Results_integration_', filename, '.mat'],...
        'corr_d_avg', 'corr_d_se', 'corr_d_fitpar');
    subplot(N_lA, N_lC, (i - 1) * N_lC + j); hold on;
    plot(d_BinCenter' * d_unit, corr_d_avg(:, end),...
        'color', 'b', 'linewidth', 1, 'marker', '.', 'markersize', 5);
    par = corr_d_fitpar(:, end, 1); par(2) = par(2) * d_unit; sigma_corr(i, j) = par(2);
    plot(d_fit, Exp2(d_fit, par(1), par(2), par(3)),...
        'linestyle', '--', 'color', [0.5 0 1], 'linewidth', 0.5);
    plot([0 100], [0 0], 'k--');
    axis([0 100 -0.1 1]); set(gca, 'XTick', 0: 10: 100, 'YTick', [-0.1 0: 0.25: 1],...
        'XTickLabel', {'0', '', '20', '', '40', '', '60', '', '80', '', '100'},...
        'YTickLabel', {'', '0', '0.25', '0.5', '0.75', '1'}); grid on;
    if j ~= 1, set(gca, 'YTickLabel', []); end
    if i ~= length(lC_list), set(gca, 'XTickLabel', []); end
    ax = gca; ax.FontSize = 6;    % for ticks
    %
    text(95, 0.85, ['\sigma = ', num2str(par(2), '%.2f')],...
        'HorizontalAlignment', 'right', 'FontSize', 7);
    text(95, 0.55, '(\mum)', 'HorizontalAlignment', 'right', 'FontSize', 7);
    if (i == N_lC) & (j == 1)
        xlabel({'Cortical distance (\mum)',...
            ['(', num2str(d_unit), ' \mum per grid)']}, 'FontSize', 8);
    end
    if j == 1, ylabel(['lA = ', num2str(lA)], 'FontSize', 10); end
    if i == 1, title(['lC = ', num2str(lC)], 'fontweight', 'normal', 'FontSize', 10); end
end
end
%
pause(2);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lC_vs_lA_supp.png']);
close;


lC_list_fine = 0.8: 0.05: 1.8;
lA_list_fine = 1.75: -0.125: 0.75;
[x, y] = meshgrid(lC_list, lA_list);
[x_fine, y_fine] = meshgrid(lC_list_fine, lA_list_fine);
sigma_corr_fine = interp2(x, y, sigma_corr, x_fine, y_fine);
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 0.8]);
subplot(1, 2, 1); imagesc(sigma_corr);
set(gca, 'XTick', 1: N_lC, 'XTickLabel', lC_list,...
    'YTick', 1: N_lA, 'YTickLabel', lA_list);
subplot(1, 2, 2); imagesc(sigma_corr_fine);
idxC = find(ismember(lC_list_fine, lC_list));
idxA = find(ismember(lA_list_fine, lA_list));
set(gca, 'XTick', idxC, 'XTickLabel', lC_list_fine(idxC),...
    'YTick', idxA, 'YTickLabel', lA_list_fine(idxA));
for k = 1: 2
    subplot(1, 2, k);
    axis square; h = colorbar;
    clim([5 11]); % clim([floor(min(sigma_corr(:))) ceil(max(sigma_corr(:)))]);
    xlabel('lC'); ylabel('lA'); title(h, '\sigma (\mum)');
end
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_N_sqrt_', num2str(N_sqrt), '_lC_vs_lA.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/ParScan_lC_vs_lA.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;


%% THEORY
dir_0 = 'D:/Study/CompNeuro/Projects/L4_Oja';
dir_2D = [dir_0, '/sim_theory_2D_fin'];
addpath(dir_0); addpath(dir_2D);
addpath(genpath('D:/Study/CompNeuro/Projects/Functions_simul'));

N_sqrt = 21; N = N_sqrt ^ 2;
rho_deg = 0.5;
sigma_A_input_deg_0 = 3.5;

x_sqrt_fine = -(N_sqrt - 1) / 2: 0.1: (N_sqrt - 1) / 2;    % unit: grid
N_sqrt_fine = length(x_sqrt_fine);
idx_xcorr2_ctr = (N_sqrt_fine - 1) / 2 + (1: N_sqrt_fine);
[jj, ii] = meshgrid(x_sqrt_fine);
r_grid_fine = sqrt(jj .^ 2 + ii .^ 2); clear ii jj
idx_ctr_fine = (N_sqrt_fine + 1) / 2;
r_fine = diag(r_grid_fine); r_fine = r_fine(idx_ctr_fine: end)';
rho_fine = r_fine / N_sqrt; drho_fine = rho_fine(2) - rho_fine(1);
idx_rho_fine_kept = find(rho_fine <= 10 / N_sqrt);
%
filename0 = ['N_sqrt_', num2str(N_sqrt), '_rhodeg', num2str(rho_deg),...
    '_sA', num2str(sigma_A_input_deg_0), '_aC1_aK1'];
load([dir_2D, '/results/Results_', filename0, '.mat'],...
    'rho_mag', 'Gn_double', 'A_rec', 'kappa');
d_unit = (rho_mag * rho_deg);
n_fine = (0: 0.0001: N_sqrt) / N_sqrt;
s_ring = 0.05;
%
Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;
IniVal = [1, 10, 0]; par_lb = [1e-5, 1e-5, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none');
% x_fit = 0: 0.1: 100;
%
lC_list_fine = 1.8: -0.02: 0.8;
lK_list_fine = 0.8: 0.02: 1.8;
ip_d_fitpar_tot = NaN(length(lC_list_fine), length(lK_list_fine), 2, 3);

tic;
for i = 1: length(lC_list_fine)
for j = 1: length(lK_list_fine)
    lC = lC_list_fine(i);
    lK = lK_list_fine(j);
    %
    par_DoG_C = [[1.493, 0.495], lC * [0.893, 1.587] / rho_deg];
    C1 = par_DoG_C(1); C2 = par_DoG_C(2);
    sigma_C_E = par_DoG_C(3); sigma_C_I = par_DoG_C(4);
    rho_C = sqrt( log((C1 * sigma_C_E ^ 4) / (C2 * sigma_C_I ^ 4)) /...
        (2 * (pi^2) * (sigma_C_E ^ 2 - sigma_C_I ^ 2)) );
    %
    sigma_grid = [[5 25] * lK; 145, 110]' / d_unit;  % row E/I, col N/B
    Wn = A_rec(1) * Gn_double(n_fine, kappa(1), sigma_grid(1, 1), sigma_grid(1, 2))...
        + A_rec(2) * Gn_double(n_fine, kappa(2), sigma_grid(2, 1), sigma_grid(2, 2));
    Kn = 1 ./ (1 - Wn); [Kn_max, idx] = max(Kn); rho_K = n_fine(idx); clear idx
    %
    rho_1 = (rho_K - rho_C) * N_sqrt;
    rho_2 = (rho_K + rho_C) * N_sqrt;
    % double-ring function
    b_func = exp(-(r_grid_fine - rho_1).^2 / (2 * s_ring^2)) +...
        exp(-(r_grid_fine - rho_2).^2 / (2 * s_ring^2));
    b_func_crr = xcorr2(b_func, b_func);
    b_func_crr = b_func_crr(idx_xcorr2_ctr, idx_xcorr2_ctr);
    b_func_crr_sqr = b_func_crr .^ 2;
    b_func_crr_1D = diag(b_func_crr_sqr); b_func_crr_1D = b_func_crr_1D(idx_ctr_fine: end)';
    % inverse Hankel transform
    b_func_crr_sqr_iHk = sum(2 * pi * drho_fine * bsxfun(@times,...
        rho_fine(idx_rho_fine_kept)' .* b_func_crr_1D(idx_rho_fine_kept)',...
        besselj(0, 2 * pi * rho_fine(idx_rho_fine_kept)' * r_fine)), 1);
    %
    ip_d_fitpar = NaN(2, 3);
    Err = @(par) Exp2(r_fine * d_unit, par(1), par(2), par(3)) - b_func_crr_sqr_iHk;
    [par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
    ip_d_fitpar(1, :) = par;
    CI = nlparci(par, residual, 'jacobian', Jacobian);
    ip_d_fitpar(2, :) = (CI(:, 2) - CI(:, 1))' / 2;
    ip_d_fitpar_tot(i, j, :, :) = reshape(ip_d_fitpar, [1 1 2 3]);
    %
    clear lC lK par_DoG_C C1 C2 sigma_C_E sigma_C_I rho_C sigma_grid Wn Kn
    clear Kn_max rho_K rho_1 rho_2 b_func b_func_crr b_func_crr_sqr b_func_crr_1D 
    clear b_func_crr_sqr_iHk ip_d_fitpar Err par residual Jacobian CI
end
%
fprintf(['i = ', num2str(i), ' / ', num2str(length(lC_list_fine)), ', ',...
    num2str(toc / 60, '%.2f'), ' min.\n']);
end
clear i j Exp2 IniVal options rho_mag rho_deg Gn_double A_rec kappa n_fine s_ring
clear x_sqrt_fine N_sqrt_fine idx_xcorr2_ctr r_grid_fine idx_ctr_fine
clear r_fine rho_fine drho_fine idx_rho_fine_kept

% save([dir_2D, '/results/ParScan_N_sqrt_', num2str(N_sqrt), '/Theory_lC_lK.mat']);
load([dir_2D, '/results/ParScan_N_sqrt_', num2str(N_sqrt), '/Theory_lC_lK.mat']);

sigma_corr_theo = ip_d_fitpar_tot(:, :, 1, 2);
%
lK_list_ff = 0.8: 0.005: 1.8;
lC_list_ff = 1.8: -0.005: 0.8;
[x, y] = meshgrid(lK_list_fine, lC_list_fine);
[x_ff, y_ff] = meshgrid(lK_list_ff, lC_list_ff);
sigma_corr_theo_ff = interp2(x, y, sigma_corr_theo, x_ff, y_ff);
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 1, 0.8]); 
subplot(1, 2, 1); imagesc(sigma_corr_theo);
set(gca, 'XTick', 1: 5: length(lK_list_fine),...
    'XTickLabel', lK_list_fine(1: 5: length(lK_list_fine)),...
    'YTick', 1: 5: length(lC_list_fine),...
    'YTickLabel', lC_list_fine(1: 5: length(lC_list_fine)));
subplot(1, 2, 2); imagesc(sigma_corr_theo_ff);
set(gca, 'XTick', 1: 20: length(lK_list_ff),...
    'XTickLabel', lK_list_ff(1: 20: length(lK_list_ff)),...
    'YTick', 1: 20: length(lC_list_ff),...
    'YTickLabel', lC_list_ff(1: 20: length(lC_list_ff)));
for k = 1: 2
    subplot(1, 2, k);
    axis square; h = colorbar;
    clim([5 13]); % clim([floor(min(sigma_corr(:))) ceil(max(sigma_corr(:)))]);
    xlabel('lK'); ylabel('lC'); title(h, '\sigma (\mum)');
end
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_Theory_lK_vs_lC.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/ParScan_Theory_lK_vs_lC.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;



figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.03, 0, 0.35, 1]); 
for k = 1: length(lC_list)
    subplot(length(lC_list), 1, k); hold on;
    scatter(lK_list, sigma_corr_sim(...
        abs(lC_list - lC_list(k)) < 1e-10, :), 35, 'r', 'filled');
    plot(lK_list_fine, sigma_corr_theo(...
        abs(lC_list_fine - lC_list(k)) < 1e-10, :), 'color', 'b', 'linewidth', 1.5);
    ylabel(['lC = ', num2str(lC_list(k))]); ylim([4 14]); set(gca, 'YTick', 5: 2: 13);
end
subplot(length(lC_list), 1, length(lC_list)); xlabel('lK');
%
pause(1);
print(gcf, '-dpng', [dir_2D, '/figures/', filename0,...
    '/ParScan_Theory_Sim_lK_vs_lC.png']);
exportgraphics(gcf, [dir_0, '/paper_writting/Figures/ai/ParScan_Theory_Sim_lK_vs_lC.eps'],...
    'BackgroundColor', 'none', 'ContentType', 'vector');
close;
