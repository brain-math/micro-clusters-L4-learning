function [Wt_BslJ_E_cpx_coeff, Corr_FT, Corr_d_iHk] =...
    func_corr_d_theory(Wt_ctr, sigma_R, rho_C, m_list, r1_iHk)
% VERY FAST

N = size(Wt_ctr, 1); N_sqrt = sqrt(N);
x_sqrt = -(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2;    % unit: grid
[j_grid, i_grid] = meshgrid(x_sqrt);
r_grid = sqrt(j_grid .^ 2 + i_grid .^ 2);
k_grid = r_grid / N_sqrt;
theta_grid = atan2(i_grid, j_grid);
idx_ctr = (N_sqrt + 1) / 2;

% Bases along alpha: Envelope & Circular Harmonics
Gau_Env_func = exp(- r_grid .^ 2 / (2 * sigma_R ^ 2));
BslJ_E_cpx = @(m)...
    besselj(abs(m), 2 * pi * rho_C * r_grid) .* exp(sqrt(-1) * m * theta_grid) .* Gau_Env_func;
% m_list = -4: 4;
Nm_list = length(m_list);    % 9
BslJ_E_base_cpx = NaN(N, Nm_list);
for m = 1: Nm_list, BslJ_E_base_cpx(:, m) = reshape(BslJ_E_cpx(m_list(m)), [N 1]); end; clear m
pinv_BslJ_E_base_cpx = pinv(BslJ_E_base_cpx);

% Decompose to bij
Wt_BslJ_E_cpx_coeff = (pinv_BslJ_E_base_cpx * Wt_ctr').';
% Wt_ctr_recover = real(BslJ_E_base_cpx * Wt_BslJ_E_cpx_coeff.')';
% RF_2D_visualize_fast(Wt_ctr, 1, -0.015, 0.015, Clr_RWB);
% RF_2D_visualize_fast(Wt_ctr_recover, 1, -0.015, 0.015, Clr_RWB);
%
% This is EXACTLY bij (complex).
Wt_BslJ_E_cpx_coeff_ft2 = NaN(N_sqrt, N_sqrt, Nm_list);
for m = 1: Nm_list
    c_tmp = reshape(Wt_BslJ_E_cpx_coeff(:, m), [N_sqrt N_sqrt]);
    Wt_BslJ_E_cpx_coeff_ft2(:, :, m) =...
        fftshift(fftn(ifftshift(c_tmp)));% / numel(c_tmp);
end
clear m c_tmp
%
% Wt_BslJ_E_cpx_coeff_2 = zeros(441, 9);
% load(savename, 'i_grid', 'j_grid');
% for i = 1: N_sqrt
% for j = 1: N_sqrt
%     tmp = exp(2 * pi * sqrt(-1) * ((i - idx_ctr) * i_grid + (j - idx_ctr) * j_grid) / N_sqrt);
%     Wt_BslJ_E_cpx_coeff_2 = Wt_BslJ_E_cpx_coeff_2 +...
%         tmp(:) * squeeze(Wt_BslJ_E_cpx_coeff_ft2(i, j, :)).';
% end
% end
% tmp = Wt_BslJ_E_cpx_coeff - Wt_BslJ_E_cpx_coeff_2; sum(tmp(:)' * tmp(:));
% clear i j i_grid j_grid tmp


% Preparation of C0 term
m_abs_list = unique(abs(m_list));    % 0: 4
Nm_abs_list = length(m_abs_list);
%
% ds = 0.1; s = 0: ds: 10;
rho = (0: 0.1: N_sqrt) / N_sqrt; drho = rho(2) - rho(1);
%
% n0 = 2 * pi * rho_C;
% cst = ( 2 * pi * exp(-2 * sigma_R^2 * n0^2) ) /...
%     (4 * pi^2 * sigma_R^4 * n0^2);
C0 = NaN(Nm_abs_list * ones(1, 4));
for i1 = 1: Nm_abs_list
m_abs_1 = m_abs_list(i1);
for i2 = 1: Nm_abs_list
m_abs_2 = m_abs_list(i2);
for i3 = 1: Nm_abs_list
m_abs_3 = m_abs_list(i3);
for i4 = 1: Nm_abs_list
m_abs_4 = m_abs_list(i4);
%
% tmp = s .* exp((-2 * s.^2) / (sigma_R^2 * n0^2)) .*...
    % besseli(m_abs_1, s) .* besseli(m_abs_2, s) .*...
    % besseli(m_abs_3, s) .* besseli(m_abs_4, s);
% C0(i1, i2, i3, i4) = ds * sum(tmp) *...
    % (sqrt(-1) ^ (- m_abs_1 + m_abs_2 + m_abs_3 - m_abs_4));
itg = rho .* exp(-8 * pi^2 * sigma_R^2 * (rho_C^2 + rho.^2)) .*...
    besseli(m_abs_1, 4 * pi^2 * sigma_R^2 * rho_C * rho) .*...
    besseli(m_abs_2, 4 * pi^2 * sigma_R^2 * rho_C * rho) .*...
    besseli(m_abs_3, 4 * pi^2 * sigma_R^2 * rho_C * rho) .*...
    besseli(m_abs_4, 4 * pi^2 * sigma_R^2 * rho_C * rho);
C0(i1, i2, i3, i4) = 4 * pi^2 * sigma_R^2 *...
    (sqrt(-1) ^ (- m_abs_1 + m_abs_2 + m_abs_3 - m_abs_4)) *...
    sum(itg) * drho;
end
end
end
end
% C0 = C0 * cst;
clear n0 cst
clear m_abs_list Nm_abs_list ds s drho rho
clear i1 i2 i3 i4 m_abs_1 m_abs_2 m_abs_3 m_abs_4 tmp


% Go loop over m
delta_m_list = m_list' - m_list;
[~, ~, m0_table] = unique(delta_m_list);
N_m0 = 2 * Nm_list - 1;
% all 'm' later is index from 1 to 9.
%
idx_xcorr2_ctr = (N_sqrt - 1) / 2 + (1: N_sqrt);
%
Corr_FT = zeros(N_sqrt, N_sqrt);

for m0_idx = 1: N_m0
idx = find(m0_table == m0_idx); n_sqrt_list = length(idx);
m1_list = mod(idx - 1, Nm_list) + 1; m2_list = ceil(idx / Nm_list); clear idx
% indices are [1;2;3;4;5] and [5;6;7;8;9]
% m values are [-4;-3;-2;-1;0] and [0;1;2;3;4] -- so m1 - m2 / m3 - m4 = m0
%
for m1_idx = 1: n_sqrt_list
m1 = m1_list(m1_idx); m2 = m2_list(m1_idx);
for m3_idx = 1: n_sqrt_list
m3 = m1_list(m3_idx); m4 = m2_list(m3_idx);
%
tmp1 = conj(xcorr2(Wt_BslJ_E_cpx_coeff_ft2(:, :, m1),...
    Wt_BslJ_E_cpx_coeff_ft2(:, :, m2)));
tmp1 = tmp1(idx_xcorr2_ctr, idx_xcorr2_ctr);
%
tmp2 = conj(xcorr2(Wt_BslJ_E_cpx_coeff_ft2(:, :, m3),...
    Wt_BslJ_E_cpx_coeff_ft2(:, :, m4)));
tmp2 = tmp2(idx_xcorr2_ctr, idx_xcorr2_ctr);
%
tmp3 = C0(abs(m_list(m1)) + 1, abs(m_list(m2)) + 1,...
    abs(m_list(m3)) + 1, abs(m_list(m4)) + 1);
Corr_FT = Corr_FT + tmp3 * (tmp1 .* tmp2);
end
end
end
clear m1_idx m3_idx m1 m2 m3 m4 tmp1 tmp2 tmp3 m0_idx n_sqrt_list m1_list m2_list


% % First look at how 2D abs(Corr_FT) over k
% k1_grid = unique(r_grid(:)) / N_sqrt;
% k1_grid_edge = [k1_grid(1); k1_grid(1: end - 1) + diff(k1_grid) / 2; k1_grid(end)];
% [Ck_abs_1_avg, Ck_abs_1_se, ~] = histogram_mean_sem(...
%     abs(Corr_FT_0(:)), r_grid(:) / N_sqrt, k1_grid_edge);
% %
% figure;
% errorbar(k1_grid(2: end) * N_sqrt, Ck_abs_1_avg(2: end), Ck_abs_1_se(2: end));
% xlim([0 10]);


% Though Corr_FT(k1, k2) has real and imag parts,
  % direct fftshift(ifft2(ifftshift( ))) turns out to be bullshit...
% Incorrect but 'good-looking' way: iHankel(?) on abs(?) ?!!
%
% r1_iHk = (0: 0.1: 10)';
% Corr_d_iHk = sum(bsxfun(@times, k_grid(:) .* abs(Corr_FT_0(:)),...
    % besselj(0, 2 * pi * k_grid(:) * r1_iHk')), 1)';
dk = k_grid(idx_ctr, idx_ctr + 1) - k_grid(idx_ctr, idx_ctr);
Corr_d_iHk = 2 * pi * (dk ^ 2) *...
    sum(bsxfun(@times, k_grid(:) .* abs(Corr_FT(:)),...
    besselj(0, 2 * pi * k_grid(:) * r1_iHk')), 1)';
Corr_d_iHk = Corr_d_iHk + real(Corr_FT(idx_ctr, idx_ctr)) / N;
%
% figure; plot(r1 * d_unit, Corr_d_iHk);
