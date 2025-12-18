function [corr_d_avg, corr_d_se, corr_d_fitpar, ip_d_avg, ip_d_se, ip_d_fitpar] =...
    func_2D_W_stat_simple_CorrMultip(r_ffwd_sqr_rounds, dist_circ, d_BinCenter, idx_fit_valid)

N_loop = size(r_ffwd_sqr_rounds, 3);

N_BinCenter = length(d_BinCenter);
d_BinEdge = [d_BinCenter(1) - (d_BinCenter(2) - d_BinCenter(1)) / 2,...
    mean([d_BinCenter(1: N_BinCenter - 1); d_BinCenter(2: N_BinCenter)], 1),...
    d_BinCenter(end) + (d_BinCenter(end) - d_BinCenter(end - 1)) / 2];
%
dist_grid = triu_new(dist_circ, 0, 1);
idx_valid = find(dist_grid < d_BinEdge(end));
dist_grid_N = repmat(dist_grid(idx_valid), [1 N_loop]);

corr_ffwd_sqr = NaN(size(dist_grid_N));
ip_ffwd_sqr = NaN(size(dist_grid_N));
for loop_k = 1: N_loop
    corr_ffwd_sqr_k = triu_new(...
        corrcoef(r_ffwd_sqr_rounds(:, :, loop_k)'), 0, 1);
    corr_ffwd_sqr(:, loop_k) = corr_ffwd_sqr_k(idx_valid);
    %
    ip_ffwd_sqr_k = triu_new(...
        r_ffwd_sqr_rounds(:, :, loop_k) * r_ffwd_sqr_rounds(:, :, loop_k)', 0, 1);
    ip_ffwd_sqr(:, loop_k) = ip_ffwd_sqr_k(idx_valid);
end

corr_d_avg = NaN(N_BinCenter, N_loop + 1);    % first 5 indiv, last one integ
corr_d_se = NaN(N_BinCenter, N_loop + 1);
ip_d_avg = NaN(N_BinCenter, N_loop + 1);    % first 5 indiv, last one integ
ip_d_se = NaN(N_BinCenter, N_loop + 1);
for loop_k = 1: N_loop
    [corr_d_avg(:, loop_k), corr_d_se(:, loop_k), ~] =...
        histogram_mean_sem(corr_ffwd_sqr(:, loop_k), dist_grid_N(:, 1), d_BinEdge);
    [ip_d_avg(:, loop_k), ip_d_se(:, loop_k), ~] =...
        histogram_mean_sem(ip_ffwd_sqr(:, loop_k), dist_grid_N(:, 1), d_BinEdge);
end
[corr_d_avg(:, end), corr_d_se(:, end), ~] =...
    histogram_mean_sem(corr_ffwd_sqr(:), dist_grid_N(:), d_BinEdge);
[ip_d_avg(:, end), ip_d_se(:, end), ~] =...
    histogram_mean_sem(ip_ffwd_sqr(:), dist_grid_N(:), d_BinEdge);

Exp2 = @(x, A, sigma, b) A * exp(- x .^ 2 / (2 * (sigma ^ 2))) + b;
IniVal = [1, 10, 0]; par_lb = [1e-5, 1e-5, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none');
%
corr_d_fitpar = NaN(3, N_loop + 1, 2);    % mean / err
ip_d_fitpar = NaN(3, N_loop + 1, 2);
for loop_k = 1: N_loop + 1
    Err_corr_k = @(par) (Exp2(d_BinCenter(idx_fit_valid)', par(1), par(2), par(3)) -...
        corr_d_avg(idx_fit_valid, loop_k)) ./ corr_d_se(idx_fit_valid, loop_k);
    [par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err_corr_k, IniVal, par_lb, par_ub, options);
    corr_d_fitpar(:, loop_k, 1) = par';
    CI = nlparci(par, residual, 'jacobian', Jacobian);
    corr_d_fitpar(:, loop_k, 2) = (CI(:, 2) - CI(:, 1)) / 2;
    %
    Err_ip_k = @(par) (Exp2(d_BinCenter(idx_fit_valid)', par(1), par(2), par(3)) -...
        ip_d_avg(idx_fit_valid, loop_k)) ./ ip_d_se(idx_fit_valid, loop_k);
    [par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err_ip_k, IniVal, par_lb, par_ub, options);
    ip_d_fitpar(:, loop_k, 1) = par';
    CI = nlparci(par, residual, 'jacobian', Jacobian);
    ip_d_fitpar(:, loop_k, 2) = (CI(:, 2) - CI(:, 1)) / 2;
end
