function [q_fin, r_fin, Eigval_fin, R2_fin, Eigval_itr, R2_itr] =...
    func_eigvec_FinelyTune(r0, Q_history, A_mat_in, K_mat, C_mat, N_itr_max, eigval_itr_thr)
% (return_val)  (eigval_tgt, eigval_tgterr_thr)
% 3 matrices: (N, N); r0: (N, N) or (N^2, 1)
N = size(A_mat_in, 1);
r0 = r0(:);    % (N^2, 1)

% ACK_mat = diag(A_mat_in(:)) * kron(C_mat, K_mat);
% ACK_mat_sym = sqrt(diag(A_mat_in(:))) * kron(C_mat, K_mat) * sqrt(diag(A_mat_in(:)));
% Operator_ACK = @(W) A_mat_in .* (K_mat * W * C_mat);
Operator_ACK_sym = @(W) sqrt(A_mat_in) .* (K_mat * (sqrt(A_mat_in) .* W) * C_mat);
if isempty(Q_history)
    Operator_ACK_sym_vec = @(w)...
        reshape(Operator_ACK_sym(reshape(w, [N N])), [N^2, 1]);
else
    Operator_ACK_sym_vec = @(w)...
        reshape(Operator_ACK_sym(reshape(w - Q_history * (Q_history' * w), [N N])), [N^2, 1]);
end

L2_norm = @(w) w / sqrt(sum(w(:) .^ 2));

q0 = L2_norm((A_mat_in(:) .^ (-1/2)) .* r0);

Eigval_itr = NaN(N_itr_max + 1, 1);
% Inner_product_0_itr = NaN(N_itr + 1, 1);
R2_itr = NaN(N_itr_max + 1, 1);
Eigval_itr(1) = sqrt(sum(Operator_ACK_sym_vec(q0) .^ 2));
% Inner_product_0_itr(1) = 1;
R2_itr(1) = 1;

q = q0; q0_avg = mean(q0);
% return_val = 0;    % failure by default.
for itr_k = 2: N_itr_max + 1
    Aq = Operator_ACK_sym_vec(q);
    q = L2_norm(Aq);
    Eigval_itr(itr_k) = sqrt(sum(Aq .^ 2));
    % Inner_product_0_itr(itr_k) = q' * q0;
    R2_itr(itr_k) = 1 - sum((q - q0) .^ 2) / sum((q0 - q0_avg) .^ 2);
    if (itr_k > 5) & (abs(Eigval_itr(itr_k - 1) - Eigval_itr(itr_k)) <= eigval_itr_thr)
        break;
    end
    %if (abs(Eigval_itr(itr_k) - eigval_tgt) <= eigval_tgterr_thr)
    %    return_val = 1;    % converge
    %    break;
    %elseif sign(Eigval_itr(itr_k - 1) - eigval_tgt) * sign(Eigval_itr(itr_k) - eigval_tgt) == -1
    %    return_val = 2;    % cross -- not really stable
    %    break;
    %end
end
q_fin = q;
r_fin = L2_norm((A_mat_in(:) .^ (1/2)) .* q);
Eigval_fin = Eigval_itr(itr_k);
% Inner_product_0_fin = Inner_product_0_itr(itr_k);
R2_fin = R2_itr(itr_k);
if itr_k ~= N_itr_max + 1
    Eigval_itr(itr_k + 1: end) = Eigval_fin;
    R2_itr(itr_k + 1: end) = R2_fin;
end