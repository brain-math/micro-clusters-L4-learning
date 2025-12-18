function RF_2D_visualize_fast(RF, d_boundary, clim_min, clim_max, Clrmap)
% RF: (Nx, Na)

[Nx, Na] = size(RF);
Nx_sqrt = sqrt(Nx); Na_sqrt = sqrt(Na);
N2 = Nx_sqrt * Na_sqrt + (Nx_sqrt - 1) * d_boundary;
% N2 = N + (N_sqrt - 1) * d_boundary;

RF4 = reshape(RF, [Nx_sqrt, Nx_sqrt, Na_sqrt, Na_sqrt]);
RF_vis = NaN(N2, N2);
for i = 1: Nx_sqrt
for j = 1: Nx_sqrt
    ifill = (i - 1) * (Na_sqrt + d_boundary) + [1: Na_sqrt];
    jfill = (j - 1) * (Na_sqrt + d_boundary) + [1: Na_sqrt];
    RF_vis(ifill, jfill) = squeeze(RF4(i, j, :, :));
end
end

figure('position', [75, 0, 900, 1000]);
axes('Units', 'Normalize', 'Position', [0 0 1 1]);    % No boundary
h = imagesc(RF_vis); set(h, 'AlphaData', ~isnan(RF_vis));    % NaN transparent
axis square;
if ~isequal(Clrmap, 0), colormap(Clrmap); end; clim([clim_min clim_max]);
axis off;