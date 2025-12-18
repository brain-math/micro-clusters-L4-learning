function [nC_max_fine, r_ffwd_sqr] =...
    func_2D_W_stat_simple_rsp(Wt, DoG_func, par_DoG_C, dtheta_deg, dphi_deg)

N = size(Wt, 1); N_sqrt = sqrt(N);

% coordinate
x_sqrt = -(N_sqrt - 1) / 2: 1: (N_sqrt - 1) / 2;    % unit: grid
[j_grid, i_grid] = meshgrid(x_sqrt);

% generate grating
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
Wt_ctr = circmat_align_2D(Wt, 1);
%
r_ffwd_sqr = NaN(N, N_theta);
for theta_i = 1: N_theta
    stim_grating_ij = stim_grating(:, :, theta_i);    % (N, N_phi);
    rL = Wt_ctr * stim_grating_ij;
    rLp_sqr = (rL .* (rL > 0)) .^ 2;
    r_ffwd_sqr(:, theta_i) = max(rLp_sqr, [], 2);
end
