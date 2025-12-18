function [nC_max_fine, r_ffwd_sqr, corr_d_ffwd_sqr_avg, corr_d_ffwd_sqr_se,...
    corr_d_fitpar_ffwd_sqr, theta_pref_ffwd_sqr, gOSI_ffwd_sqr,...
    r_rec, corr_d_rec_avg, corr_d_rec_se,...
    corr_d_fitpar_rec, theta_pref_rec, gOSI_rec] =...
    func_2D_W_stat_rho0(W_trace_select_k,...
    K_mat, DoG_func, par_DoG_C, dtheta_deg, dphi_deg, dist_grid, d_BinCenter, idx_fit_valid)
% W_trace_select_k = W_trace(:, :, trace_idx_select);    % (N, N, N_t)
% dist_grid = dist_circ;

% if idx_fit_valid = [], then do not fit

N = size(W_trace_select_k, 1); N_sqrt = sqrt(N);
N_t = size(W_trace_select_k, 3);

% coordinate
x_sqrt = -(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2;    % unit: grid
[j_grid, i_grid] = meshgrid(x_sqrt);

% generate grating
% See supp_math_absFT_as_grating_rsp.jpg -- absFT of W is the linear ffwd grating response
  % (of each (nx, ny), i.e. (n cos(theta), n sin(theta) ) -- max over all phase.
  % so W_absFT_trace is 'tuning curve' -- 2nd dim as stim (nx, ny).
  % And W_absFT_trace^2 is linear response then go through ()_+^2 -- easier for theory.
% We can use this for analytical analysis, but don't use it for numerical Corr(d)
  % and pref ori map -- n too rough, esp when N_sqrt is small.
% So to do these numerically we'll use real grating.
%
% dtheta_deg = 15;
theta_stim = [dtheta_deg: dtheta_deg: 180] * (pi / 180); N_theta = length(theta_stim);
% dphi_deg = 15;
phi_stim = [dphi_deg: dphi_deg: 360] * (pi / 180); N_phi = length(phi_stim);
% For n: Let's use single, optimal n.
dr_fine = 0.01; r_fine = [0: dr_fine: (N_sqrt + 1) / 2];
nr_fine = r_fine / N_sqrt;
C_func_r_fine = DoG_func(r_fine, par_DoG_C);
itg = bsxfun(@times, (2 * pi * dr_fine) * r_fine .* C_func_r_fine,...
    besselj(0, 2 * pi * (nr_fine' * r_fine)));
C_func_nr_fine = sum(itg, 2)';    % Hankel transformation
[~, idx] = max(C_func_nr_fine); nC_max_fine = nr_fine(idx);    % * N_sqrt
clear dr_fine r_fine nr_fine C_func_r_fine itg idx C_func_nr_fine
%
% figure; k = 1;
% for theta_i = 1: N_theta
%     theta = theta_list(theta_i);
%     grating_2_test = cos(2 * pi * nC_max_fine * (x_grid * cos(theta) + y_grid * sin(theta)));
%     subplot(1, N_theta, k);
%     pcolor(x_grid, y_grid, grating_2_test, 'edgecolor','none');
%     axis square; clim([-1 1]);
%     k = k + 1;
% end
%
stim_grating = NaN(N, N_phi, N_theta);
for theta_i = 1: N_theta
    theta = theta_stim(theta_i);
    for phi_k = 1: N_phi 
        stim_grating(:, phi_k, theta_i) = cos(2 * pi * nC_max_fine *...
            (j_grid(:) * cos(theta) + i_grid(:) * sin(theta)) + phi_stim(phi_k));
    end
end


% neuronal response
% circshift W to center will remove the non-zero Corr(d) asympt value numerically.
% Don't ask me why but it's necessary.
idx_ctr = (N_sqrt + 1) / 2;
W_trace_select_ctr = NaN(N, N, N_t);
for ti = 1: N_t
for xi = 1: N
    i = mod(xi - 1, N_sqrt) + 1; j = ceil(xi / N_sqrt);
    Wk = reshape(W_trace_select_k(xi, :, ti), [N_sqrt, N_sqrt]);
    W_trace_select_ctr(xi, :, ti) = reshape(circshift(Wk, [idx_ctr - i, idx_ctr - j]), [1 N]);
end
end
%
% r_ffwd_linear = NaN(N, N_theta, N_t);
r_ffwd_sqr = NaN(N, N_theta, N_t);
r_rec = NaN(N, N_theta, N_t);
for theta_i = 1: N_theta
    stim_grating_ij = stim_grating(:, :, theta_i);    % (N, N_phi);
    for t_k = 1: N_t
        rL = W_trace_select_ctr(:, :, t_k) * stim_grating_ij;
        rLp_sqr = (rL .* (rL > 0)) .^ 2;
        r_rec_k = K_mat * rLp_sqr;
        %
        % r_ffwd_linear(:, theta_i, t_k) = max(rL, [], 2);
        r_ffwd_sqr(:, theta_i, t_k) = max(rLp_sqr, [], 2);
        r_rec(:, theta_i, t_k) = max(r_rec_k, [], 2);
        clear rL rLp_sqr
    end
end


% thete_pref and gOSI
% r_3d: (N, N_theta, N_t)
theta_pref_func = @(r_3d, theta_stim) (1 / 2) * (pi + atan2(...
    squeeze(mean(bsxfun(@times, r_3d, sin(pi + 2 * theta_stim)), 2)),...
    squeeze(mean(bsxfun(@times, r_3d, cos(pi + 2 * theta_stim)), 2)) ));    % (N, N_t)
theta_pref_ffwd_sqr = theta_pref_func(r_ffwd_sqr, theta_stim);
theta_pref_rec = theta_pref_func(r_rec, theta_stim);
% Avg_sin = squeeze(mean(bsxfun(@times, r_ffwd_sqr, sin(pi + 2 * theta_stim)), 2));
% Avg_cos = squeeze(mean(bsxfun(@times, r_ffwd_sqr, cos(pi + 2 * theta_stim)), 2));
% theta_pref_ffwd_sqr = (atan2(Avg_sin, Avg_cos) + pi) / 2;
%
gOSI_func = @(r_3d)...
    squeeze(abs(sum(bsxfun(@times, r_3d, exp(2 * sqrt(-1) * theta_stim)), 2))) ./...
    squeeze(sum(r_3d, 2));
gOSI_ffwd_sqr = gOSI_func(r_ffwd_sqr);
gOSI_rec = gOSI_func(r_rec);
% A1 = squeeze(sum(r_ffwd_sqr, 2));
% A2 = squeeze(abs(sum(bsxfun(@times, r_ffwd_sqr, exp(2 * sqrt(-1) * theta_stim)), 2)));
% gOSI_ffwd_sqr = A2 ./ A1;

% dist and corr
% d_BinCenter = [1, sqrt(2), 2, sqrt(5), 3: 0.5: 10];
% d_BinCenter = [1, sqrt(2), 2, sqrt(5), sqrt(8), 3.5: 0.5: 10];
N_BinCenter = length(d_BinCenter);
d_BinEdge = [d_BinCenter(1) - (d_BinCenter(2) - d_BinCenter(1)) / 2,...
    mean([d_BinCenter(1: N_BinCenter - 1); d_BinCenter(2: N_BinCenter)], 1),...
    d_BinCenter(end) + (d_BinCenter(end) - d_BinCenter(end - 1)) / 2];
%
dist_grid = triu_new(dist_grid, 0, 1);    % triu_new(x, needDiag, NonNaNElementsOnly)
idx_valid = find(dist_grid < d_BinEdge(end));
%
%corr_d_ffwd_linear_avg = NaN(N_BinCenter, N_t);
%corr_d_ffwd_linear_se = NaN(N_BinCenter,  N_t);
corr_d_ffwd_sqr_avg = NaN(N_BinCenter, N_t);
corr_d_ffwd_sqr_se = NaN(N_BinCenter, N_t);
corr_d_rec_avg = NaN(N_BinCenter, N_t);
corr_d_rec_se = NaN(N_BinCenter, N_t);
%
if ~isempty(idx_fit_valid)
Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;
IniVal = [1, 10, 0]; par_lb = [1e-5, 1e-5, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none');
end
corr_d_fitpar_ffwd_sqr = NaN(3, N_t, 2);    % mean / err
corr_d_fitpar_rec = NaN(3, N_t, 2);
%
for t_k = 1: N_t
    %c_ffwd_linear = triu_new(corrcoef(r_ffwd_linear(:, :, t_k)'), 0, 1);
    %[corr_d_ffwd_linear_avg(:, t_k), corr_d_ffwd_linear_se(:, t_k), ~] =...
    %    histogram_mean_sem(c_ffwd_linear, dist_grid, d_BinEdge);
    c_ffwd_sqr = triu_new(corrcoef(r_ffwd_sqr(:, :, t_k)'), 0, 1);
    % c_ffwd_sqr = triu_new(r_ffwd_sqr(:, :, t_k) * r_ffwd_sqr(:, :, t_k)', 0, 1);
    [corr_d_ffwd_sqr_avg(:, t_k), corr_d_ffwd_sqr_se(:, t_k), ~] =...
        histogram_mean_sem(c_ffwd_sqr(idx_valid), dist_grid(idx_valid), d_BinEdge);
    if ~isempty(idx_fit_valid)
        Err1 = @(par) (Exp2(d_BinCenter(idx_fit_valid)', par(1), par(2), par(3)) -...
            corr_d_ffwd_sqr_avg(idx_fit_valid, t_k)) ./ corr_d_ffwd_sqr_se(idx_fit_valid, t_k);
        [par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err1, IniVal, par_lb, par_ub, options);
        corr_d_fitpar_ffwd_sqr(:, t_k, 1) = par';
        CI = nlparci(par, residual, 'jacobian', Jacobian);
        corr_d_fitpar_ffwd_sqr(:, t_k, 2) = (CI(:, 2) - CI(:, 1)) / 2;
    end
    clear c_ffwd_sqr Err1 par residual Jacobian CI
    %
    c_rec = triu_new(corrcoef(r_rec(:, :, t_k)'), 0, 1);
    [corr_d_rec_avg(:, t_k), corr_d_rec_se(:, t_k), ~] =...
        histogram_mean_sem(c_rec(idx_valid), dist_grid(idx_valid), d_BinEdge);
    if ~isempty(idx_fit_valid)
        Err2 = @(par) (Exp2(d_BinCenter(idx_fit_valid)', par(1), par(2), par(3)) -...
            corr_d_rec_avg(idx_fit_valid, t_k)) ./ corr_d_rec_se(idx_fit_valid, t_k);
        [par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err2, IniVal, par_lb, par_ub, options);
        corr_d_fitpar_rec(:, t_k, 1) = par';
        CI = nlparci(par, residual, 'jacobian', Jacobian);
        corr_d_fitpar_rec(:, t_k, 2) = (CI(:, 2) - CI(:, 1)) / 2;
    end
    clear c_rec Err2 par residual Jacobian CI
end
