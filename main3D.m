%% INDEXING MATRICES FOR EFFICIENT STORAGE AND EINSTEIN SUMMATION
% IND_K4 = Indexing into the symmetric parts of the elastic tangent moduli
% IND    = Indexing into the 4th-order tensor in the form (1,1,i,j), (2,1,i,j), and so on
% IND_G  = Indexing into the symmetric part of the Green operator in the Fourier space
IND_K4 = [1, 1, 1, 1;  1, 1, 1, 2;  1, 1, 1, 3;  1, 1, 2, 1;  1, 1, 2, 2;  % 1  --> 5
    1, 1, 2, 3;  1, 1, 3, 1;  1, 1, 3, 2;  1, 1, 3, 3;  1, 2, 1, 2;  % 6  --> 10
    1, 2, 1, 3;  1, 2, 2, 1;  1, 2, 2, 2;  1, 2, 2, 3;  1, 2, 3, 1;  % 11 --> 15
    1, 2, 3, 2;  1, 2, 3, 3;  1, 3, 1, 3;  1, 3, 2, 1;  1, 3, 2, 2;  % 16 --> 20
    1, 3, 2, 3;  1, 3, 3, 1;  1, 3, 3, 2;  1, 3, 3, 3;  2, 1, 2, 1;  % 21 --> 25
    2, 1, 2, 2;  2, 1, 2, 3;  2, 1, 3, 1;  2, 1, 3, 2;  2, 1, 3, 3;  % 26 --> 30
    2, 2, 2, 2;  2, 2, 2, 3;  2, 2, 3, 1;  2, 2, 3, 2;  2, 2, 3, 3;  % 31 --> 35
    2, 3, 2, 3;  2, 3, 3, 1;  2, 3, 3, 2;  2, 3, 3, 3;  3, 1, 3, 1;  % 36 --> 40
    3, 1, 3, 2;  3, 1, 3, 3;  3, 2, 3, 2;  3, 2, 3, 3;  3, 3, 3, 3]; % 41 --> 45
IND = [1,  4,  7,  2,  5,  8,  3,  6,  9;   % (1,1,i,j) -- equivalent to for j = 1:3, for i = 1:3
       4, 25, 28, 12, 26, 29, 19, 27, 30;   % (2,1,i,j)
       7, 28, 40, 15, 33, 41, 22, 37, 42;   % (3,1,i,j)
       2, 12, 15, 10, 13, 16, 11, 14, 17;   % (1,2,i,j)
       5, 26, 33, 13, 31, 34, 20, 32, 35;   % (2,2,i,j)
       8, 29, 41, 16, 34, 43, 23, 38, 44;   % (3,2,i,j)
       3, 19, 22, 11, 20, 23, 18, 21, 24;   % (1,3,i,j)
       6, 27, 37, 14, 32, 38, 21, 36, 39;   % (2,3,i,j)
       9, 30, 42, 17, 35, 44, 24, 39, 45];  % (3,3,i,j)
IND_G = [1  1; 1  2; 1  3; 2  2; 2  3; 3  3];
%% REGULAR MESH
ndim = 3;                       % problem dimension
L = ones(1,ndim);               % RVE occupies the domain [-L(1), L(1)] x [-L(2), L(2)] x [-L(3), L(3)]
n = repmat(64, [1, ndim]);      % number of nodes per each direction
h = prod(2*L./n);               % volume of one element boxing the grid point
vol = prod(2*L);                % volume of the RVE
scale = L./pi;  % for scaling the wave frequ
x1D = cell(1,ndim);     % one-dimensional coordinates in each direction
xi1D = cell(1,ndim);    % one-dimensional Fourier frequencies in each direction
for a = 1:ndim
    if rem(n(a),2) == 1
        x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2));
        xi1D{a} = [0 : fix(n(a)/2), -fix(n(a)/2) : -1] / scale(a);
    else
        x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2)-1);
        xi1D{a} = [0: fix(n(a)/2)-1, 0, -fix(n(a)/2)+1 : -1] / scale(a);
    end
end
x = zeros([n, ndim]); xi = zeros([n, ndim]);                                    % three-dimensional coordinates (x)
[x(:,:,:,1), x(:,:,:,2), x(:,:,:,3)] = ndgrid(x1D{1}, x1D{2}, x1D{3});          % and three-dimensional Fourier frequencies (xi)
[xi(:,:,:,1), xi(:,:,:,2), xi(:,:,:,3)] = ndgrid(xi1D{1}, xi1D{2}, xi1D{3});	
phase_IND = ones(n(1), n(2), n(3));                                     % matrix phase is indicated by integer 1
phase_IND( (x(:,:,:, 1) <= 0.75*L(1) & x(:,:,:, 1) >= -0.75*L(1)) ...   % rectangular inclusions are indicated by integer 2
    & ( x(:,:,:, 2) <=  0.75*L(2) & x(:,:,:, 2) >= -0.75*L(2)) ...
    & ((x(:,:,:, 3) <= -0.50*L(3) & x(:,:,:, 3) >= -0.75*L(3)) ...
    |  (x(:,:,:, 3) <=  0.75*L(3) & x(:,:,:, 3) >=  0.5*L(3))) ) = 2;
vol_frac = 0.025; radius = (6 * vol_frac * prod(L) / pi)^(1/3);         % RVE volume = 8*L(1)*L(2)*L(3)
phase_IND(sum(x.^2, ndim + 1) <= radius*radius) = 3;                    % circular inclusion is indicated by integer 3
%% MATERIAL LAWS
nu = [0.33, 0.30, 0.34];                % Poisson ratios for two material phases
E = [69, 5*69, 20*69];                  % Young moduli for two material phases
mu_ = E./(2*(1+nu));                    % Shear moduli for two material phases
lambda_ = E.*nu./((1+nu).*(1-2*nu));    % Lame parameters for two material phases
beta_ = 2*nu./(1-2*nu);                 % Material parameters for Neo-Hookean materials
constitutive_laws = cell(1,2);          % containing function handles representing various constitutive laws
constitutive_laws{1} = @(F) neo_hookean(F, mu_(1), beta_(1), IND_K4);
constitutive_laws{2} = @(F) neo_hookean(F, mu_(2), beta_(2), IND_K4);
constitutive_laws{3} = @(F) saint_venant(F, lambda_(3), mu_(3), IND_K4);
%% PRE-PROCESSING
max_iter = 20;      TOL = 1e-8;             % parameters for stopping Newton-Raphson (NR) iteration
max_iter_CG = 20;    TOL_CG = 1e-6;         % parameters for stopping the CG solver for NR each iteration
num_phases = length(constitutive_laws);     % number of phases
phase = cell(1, num_phases);                % cell of variables keeping boolean values for each material phase
for p = 1:num_phases, phase{p} = (phase_IND == p); end      % phase{p}(i,j,k) = 1 if X(i,j,k) belongs to the p-th phase
n_nodes = prod(n);                                          % total number of nodes
nDoFs_sqrt = sqrt(n_nodes*ndim*ndim);       % square root of total number of DoFs
F_Macro = [1.5, 0.25, 0; 0, 1, 0; 0, 0, 1]; % macroscopic deformation
Green = xi(:,:,:,IND_G(:,1)).*xi(:,:,:,IND_G(:,2))./sum(xi.^2, 4);	% gpuArray  % Green operator
Green(isnan(Green)) = 0;                                                        % modify elements that are divided by zeros
%% SOLVE CELL PROBLEM
F = permute(repmat(F_Macro, 1, 1, n_nodes), [3, 1, 2]);     % initialize gradient deformation
P = zeros(n_nodes,ndim,ndim);  % gpuArray                   % 1st Piola-Kirchhoff stress
K4 = zeros(n_nodes,ndim*ndim*(ndim*ndim+1)/2);  % gpuArray  % tangent moduli             
for iter = 1:max_iter
    for p = 1 : num_phases
        [P(phase{p},:,:), K4(phase{p},:)] = constitutive_laws{p}(F(phase{p},:,:));
    end
    PP = reshape(P, n(1), n(2), n(3), ndim, ndim);
    KK4 = reshape(K4, n(1), n(2), n(3), ndim*ndim*(ndim*ndim+1)/2);
    lhs = @(dF) stiff_func3D(dF, KK4, Green, IND);          % left-hand side
    rhs = residual_func3D(PP, Green);                       % right-hand side
    [dF, FLAG] = pcg(lhs, rhs, TOL_CG, max_iter);           % other possibilities: pcg, bicgstab, bicgstabl, cgs
    if dF == 0; fprintf('Step = %d -- Update vanishes\n', iter); end    
    F = F + reshape(dF, n_nodes, ndim, ndim);                   % update gradient field
    residual = norm(rhs, 2) / nDoFs_sqrt;                       % residual error
    fprintf('residual = %6.5e\n', residual);
    if residual < TOL, fprintf('Step = %d -- Residual = %6.5e\n', iter, residual); break; end
end
%% MACROSCOPIC ELASTIC TANGENT
S = zeros(n_nodes, ndim, ndim, ndim, ndim);  % gpuArray     % fluctuation sensitivities
[I, J] = ind2sub([ndim, ndim], 1:ndim*ndim);                % linear indexing -> subscript indexing 
for p = 1 : num_phases
    [P(phase{p},:,:), K4(phase{p}, :)] = constitutive_laws{p}(F(phase{p},:,:));
end
KK4 = reshape(K4, n(1), n(2), n(3), ndim*ndim*(ndim*ndim+1)/2);  % gpuArray
for i = 1:length(IND)
    lhs = @(SVect) stiff_func3D(SVect, KK4, Green, IND);
    T = reshape(KK4(:,:,:,IND(i,:)), n(1), n(2), n(3), ndim, ndim);
    rhs = residual_func3D(T, Green);
    [SVect, ~] = pcg(lhs, rhs, 1e-10, 50);
    S(:,:,:,I(i),J(i)) = reshape(SVect, n_nodes, ndim, ndim);
end
K4_Macro = zeros(1, length(IND_K4));
for i = 1:length(IND_K4)
    idx = sub2ind([ndim, ndim], IND_K4(i,3), IND_K4(i,4));  SS = squeeze(S(:,:,:,idx));
    idx = sub2ind([ndim, ndim], IND_K4(i,1), IND_K4(i,2));  
    K4_Macro(i) = h*(sum(K4(:,i)) + sum(K4(:,IND(idx,:)) .* SS(:,:), 'all')) / vol;
end
%% HELPER FUNCTIONS
function v = residual_func3D(P, G)  % Compute right-hand side for the CG solver
Phat = zeros(size(P)); ndim = 3;
for i = 1:ndim*ndim, Phat(:,:,:,i) = fftn(P(:,:,:,i)); end
v = zeros(size(P));  % gpuArray
for i = 1:ndim
    v(:,:,:,i,1) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[1, 2, 3]), 4) );
    v(:,:,:,i,2) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[2, 4, 5]), 4) );
    v(:,:,:,i,3) = -ifftn( sum(squeeze(Phat(:,:,:,i,:)).*G(:,:,:,[3, 5, 6]), 4) );
end
v = v(:);       % roll tensor into column vector
end
function v = stiff_func3D(dF_vector, K4, G, IND) % Define linear function AFUN for the CG solver
[nx, ny, nz, ~] = size(G);  ndim = 3;
dF = reshape(dF_vector, nx, ny, nz, ndim, ndim);    % unroll vector into tensor-matrix
dK = zeros(nx, ny, nz, ndim, ndim);  % gpuArray     % tangent increment
for i = 1:length(IND), dK(:,:,:,i) = fftn( sum(K4(:,:,:, IND(i,:)) .* dF(:,:,:, :), 4)); end
v = zeros(nx, ny, nz, ndim, ndim);   % gpuArray
for i = 1:ndim
    v(:,:,:,i,1) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[1, 2, 3]), 4) );
    v(:,:,:,i,2) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[2, 4, 5]), 4) );
    v(:,:,:,i,3) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[3, 5, 6]), 4) );
end
v = v(:);      % roll tensor into column vector.
end