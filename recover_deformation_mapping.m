%% RECOVER THE DEFORMATION MAPPING
F_macro_grid = permute(repmat(F_Macro, [1, 1, n]), [3, 4, 5, 1, 2]);    % macrscopic F at all grid points
F_fluc_grid = reshape(F, [n, ndim, ndim]) - F_macro_grid;               % fluctuation gradient       
xi_squared = sum(xi.^2, 4);                                             % wavenumbers squared
FF_Fluc_Fourier = zeros(size(F_fluc_grid));                             % FFT of fluctution gradient
for i = 1:ndim*ndim, FF_Fluc_Fourier(:,:,:,i) = fftn(F_fluc_grid(:,:,:,i)); end
varphi_fluc = zeros(size(x));               % fluctuation field
varphi = zeros(size(x));                    % microscopic deformation mapping
for i = 1:ndim
    What = 1i * sum(squeeze(FF_Fluc_Fourier(:,:,:,i,:)) .* xi, 4) ./ xi_squared;    What(xi_squared == 0) = 0;
    varphi_fluc(:,:,:,i) = -ifftn( What );
    varphi(:,:,:,i) = sum(squeeze(F_macro_grid(:,:,:,i,:)) .* x, 4) + varphi_fluc(:,:,:,i);
end