%% REGULAR MESH
run indexing.m; 
ndim = 3; L = ones(1,ndim); 
if ~exist('n', 'var')
    n = repmat(32, 1,3);
elseif length(n) == 1
    n = repmat(n, [1, ndim]);
end
scale = L./pi;
x1D = cell(1,ndim); xi1D = cell(1,ndim);
for a = 1:ndim
    if rem(n(a),2) == 1
        x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2));
        xi1D{a} = [0 : fix(n(a)/2), -fix(n(a)/2) : -1] / scale(a);
    else
        x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2)-1);
        xi1D{a} = [0: fix(n(a)/2)-1, 0, -fix(n(a)/2)+1 : -1] / scale(a);
    end
end
x = zeros([n, ndim]); xi = zeros([n, ndim]);
[x(:,:,:,1), x(:,:,:,2), x(:,:,:,3)] = ndgrid(x1D{1}, x1D{2}, x1D{3});
[xi(:,:,:,1), xi(:,:,:,2), xi(:,:,:,3)] = ndgrid(xi1D{1}, xi1D{2}, xi1D{3});
phase_IND = ones(n);
volfrac = 0.08;
phase_IND(sum(x.^2, ndim+1) <= (6*volfrac*prod(L)/pi)^(2/3)) = 2;

%% MATERIAL LAWS
phase_contrast = 50;
bulk_ = 0.833*[1, phase_contrast];
mu_ = 0.386*[1, phase_contrast];
lambda_ = bulk_ - 2/3 * mu_;
constitutiveLaws = cell(1,2);
constitutiveLaws{1} = @(F) saint_venant(F, lambda_(1), mu_(1), IND_K4);
constitutiveLaws{2} = @(F) saint_venant(F, lambda_(2), mu_(2), IND_K4);

%% PRE-PROCESSING
max_iter = 7;       TOL = 5e-10;
max_iter_CG = 50;   TOL_CG = 1e-6;
num_phases = length(constitutiveLaws);
phase = cell(1, num_phases);    phase_nDoFs = zeros(1, num_phases);
for p = 1:num_phases
    phase{p} = (phase_IND == p);
    phase_nDoFs(p) = length(find(phase{p}));
end
n_nodes = prod(n); nDoFs_sqrt = sqrt(n_nodes*ndim*ndim);
F_Macro = eye(3, 3);
F_Macro(1,2) = 0.4;    % Macroscopic
Green = xi(:,:,:,IND_G(:,1)).*xi(:,:,:,IND_G(:,2))./sum(xi.^2, 4);	% Green operator
Green(isnan(Green)) = 0;	% modify elements that are divided by zeros

%% SOLVE CELL PROBLEM
F = permute(repmat(F_Macro, 1, 1, n_nodes), [3, 1, 2]);  % initialize gradient deformation
P = zeros(n_nodes,ndim,ndim);                            % 1st Piola-Kirchhoff stress
T = zeros(n_nodes,ndim*ndim*(ndim*ndim+1)/2);            % tangent moduli

tic
for iter = 1:max_iter
    t = tic;
    for p = 1 : num_phases
        [P(phase{p},:,:), T(phase{p}, :)] = constitutiveLaws{p}(F(phase{p},:,:));
    end
    P_grid = reshape(P, [n, ndim, ndim]);
    K4_grid = reshape(T, [n, ndim*ndim*(ndim*ndim+1)/2]);
    lhs = @(dF) stiff_func3D(dF, K4_grid, Green, IND);        % left-hand side
    rhs = residual_func3D(P_grid, Green);               % right-hand side
    [dF, FLAG] = pcg(lhs, rhs, TOL_CG, max_iter_CG); % possibilities: pcg, cgs, bicgstab,

    F = F + reshape(dF, n_nodes, ndim, ndim);   % update gradient field
    residual = norm(rhs, 2) / nDoFs_sqrt;            % residual error
    fprintf('residual = %6.5e\n', residual);
    if residual < TOL
        fprintf('Step = %d -- Residual = %6.5e\n', iter, residual);
        break;
    end
end
toc

%% COMPUTE THE MACROSCOPIC ENERGY
E = zeros(n_nodes, ndim, ndim);
E2 = zeros(n_nodes, ndim, ndim);
for i = 1:ndim
    for j = 1:ndim
        E(:,i,j) = 0.5*(sum(F(:,:,i) .* F(:,:,j), 2) - double(i == j));
    end
end
for i = 1:ndim
    for j = 1:ndim
        E2(:,i,j) = sum(E(:,:,i) .* E(:,:,j), 2);
    end
end
traceE = E(:,1,1) + E(:,2,2) + E(:,3,3);
traceE2 = E2(:,1,1) + E2(:,2,2) + E2(:,3,3);
energy_grid_points = zeros(n_nodes,1);
for p = 1 : num_phases
    energy_grid_points(phase_IND == p) ...
        = 0.5*lambda_(p)*traceE(phase_IND == p).^2 + mu_(p)*traceE2(phase_IND == p);
end
energy_macro = 1/prod(2*L) * sum(energy_grid_points, 'all') * prod(2*L./n);

%% HELPER FUNCTIONS
% Compute right-hand side for the CG solver
function v = residual_func3D(P, G)
Phat = zeros(size(P)); ndim = 3;
for i = 1:ndim*ndim, Phat(:,:,:,i) = fftn(P(:,:,:,i)); end
v = zeros(size(P));
for i = 1:ndim
    v(:,:,:,i,1) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[1, 2, 3]), 4) );
    v(:,:,:,i,2) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[2, 4, 5]), 4) );
    v(:,:,:,i,3) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[3, 5, 6]), 4) );
end
v = v(:);       % roll tensor into column vector
end
% Define linear function AFUN for the CG solver
function v = stiff_func3D(dF_vector, K4, G, IND)
n = zeros(1,3); [n(1), n(2), n(3), ~] = size(G); ndim = 3;
dF = reshape(dF_vector, [n, ndim, ndim]);       % unroll vector into tensor-matrix
dK = zeros([n, ndim, ndim]);                    % tangent increment

for i = 1:length(IND)
    dK(:,:,:,i) = fftn( sum(K4(:,:,:, IND(i,:)) .* dF(:,:,:, :), 4)); 
end
v = zeros([n, ndim, ndim]);
for i = 1:ndim
    v(:,:,:,i,1) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[1, 2, 3]), 4) );
    v(:,:,:,i,2) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[2, 4, 5]), 4) );
    v(:,:,:,i,3) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[3, 5, 6]), 4) );
end
v = v(:);      % roll tensor into column vector.
end