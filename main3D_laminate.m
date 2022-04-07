clc%% REGULAR MESH
close all; run indexing.m
ndim = 3; L = ones(1,ndim); n = [3*3, 3, 3]; scale = L./pi;
h = prod(2*L./n); vol = prod(2*L);  x1D = cell(1,ndim); xi1D = cell(1,ndim);
for a = 1:ndim
    if rem(n(a),2) == 1
        x1D{a} = scale(a) * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2));
        xi1D{a} = [0 : fix(n(a)/2), -fix(n(a)/2) : -1] / scale(a);
    else
        x1D{a} = scale{a} * (2*pi/n(a)) * (-fix(n(a)/2) : fix(n(a)/2)-1);
        xi1D{a} = [0: fix(n(a)/2)-1, 0, -fix(n(a)/2)+1 : -1] / scale(a);
    end
end
x = zeros([n, ndim]); xi = zeros([n, ndim]);
[x(:,:,:,1), x(:,:,:,2), x(:,:,:,3)] = ndgrid(x1D{1}, x1D{2}, x1D{3});
[xi(:,:,:,1), xi(:,:,:,2), xi(:,:,:,3)] = ndgrid(xi1D{1}, xi1D{2}, xi1D{3});

phase_IND = zeros(n);
phase_IND(x(:,:,:,1) >= -    L(1) & x(:,:,:,1) < -1/3*L(1)) = 1;
phase_IND(x(:,:,:,1) >= -1/3*L(1) & x(:,:,:,1) <  1/3*L(1)) = 2;
phase_IND(x(:,:,:,1) >=  1/3*L(1) & x(:,:,:,1) <      L(1)) = 3;

%% MATERIAL LAWS
nu = [0.33, 0.30, 0.34];
E = [69, 200, 50];
mu_ = E./(2*(1+nu));
lambda_ = E.*nu./((1+nu).*(1-2*nu));
beta_ = 2*nu./(1-2*nu);
constitutiveLaws = cell(1, 3);
constitutiveLaws{1} = @(F) neo_hookean(F, mu_(1), beta_(1), IND_K4);
constitutiveLaws{2} = @(F) saint_venant(F, lambda_(2), mu_(2), IND_K4);
constitutiveLaws{3} = @(F) neo_hookean(F, mu_(3), beta_(3), IND_K4);

%% PRE-PROCESSING
max_iter = 100;      TOL = 1e-14;
max_iter_CG = 50;    TOL_CG = 1e-12;
num_phases = length(constitutiveLaws);
phase = cell(1, num_phases);
for p = 1:num_phases
    phase{p} = (phase_IND == p);
end
n_nodes = prod(n);
nDoFs_sqrt = sqrt(n_nodes*ndim*ndim);
F_Macro = [1.0, 0.4, 0.0;
    0.1, 1.0, 0.0;
    0.1, 0.0, 1.0];
Green = xi(:,:,:,IND_G(:,1)).*xi(:,:,:,IND_G(:,2))./sum(xi.^2, 4);	% Green operator
Green(isnan(Green)) = 0;                        % modify elements that are divided by zeros
%% SOLVE CELL PROBLEM
F = permute(repmat(F_Macro, 1, 1, n_nodes), [3, 1, 2]);  % initialize gradient deformation
P = zeros(n_nodes,ndim,ndim);                            % 1st Piola-Kirchhoff stress
K4 = zeros(n_nodes,ndim*ndim*(ndim*ndim+1)/2);            % tangent moduli
for iter = 1:max_iter
    for p = 1 : num_phases
        [P(phase{p},:,:), K4(phase{p}, :)] = constitutiveLaws{p}(F(phase{p},:,:));
    end
    P_grid = reshape(P, n(1), n(2), n(3), ndim, ndim);
    K4_grid = reshape(K4, n(1), n(2), n(3), ndim*ndim*(ndim*ndim+1)/2);
    lhs = @(dF) stiff_func3D(dF, K4_grid, Green, IND);        % left-hand side
    rhs = residual_func3D(P_grid, Green);               % right-hand side
    [dF, FLAG] = pcg(lhs, rhs, TOL_CG, max_iter_CG); % possibilities: pcg, bicgstab, bicgstabl, cgs
    if dF == 0
        fprintf('Step = %d -- Update vanishes\n', iter);
    end
    F = F + reshape(dF, n_nodes, ndim, ndim);    % update gradient field
    % Convergence test
    residual = norm(rhs, 2) / nDoFs_sqrt;
    if residual < TOL
        break;
    end
end

%% DISPLAY DEFORMATION GRADIENT
F = reshape(F, n(1), n(2), n(3), ndim, ndim);
format longg
fprintf('Numerical solution \n=============================\n')
disp('F[1] = ');
disp(squeeze(F(fix(n(1)/6), 1, 1, :,:)))
disp('F[2] = ');
disp(squeeze(F(fix(n(1)/6)+n(1)/3, 1, 1, :,:)))
disp('F[3] = ');
disp(squeeze(F(fix(n(1)/6)+2*n(1)/3, 1, 1, :,:)))
F = reshape(F, n_nodes, ndim, ndim);

F_analytical = readmatrix('deformationGradient.csv');
fprintf('Exact solution \n=============================\n')
disp('Fexact_{i1} = ');
disp(F_analytical(:,1))
disp('Fexact_{i2} = ');
disp(F_analytical(:,2))
disp('Fexact_{i3} = ');
disp(F_analytical(:,3))

%% MACROSCOPIC ELASTIC TANGENT
S = zeros(n_nodes, ndim, ndim, ndim, ndim);
[I, J] = ind2sub([ndim, ndim], 1:ndim*ndim);
for p = 1 : num_phases
    [P(phase{p},:,:), K4(phase{p}, :)] = constitutiveLaws{p}(F(phase{p},:,:));
end
K4_grid = reshape(K4, n(1), n(2), n(3), ndim*ndim*(ndim*ndim+1)/2);
for i = 1:length(IND)
    lhs = @(SVect) stiff_func3D(SVect, K4_grid, Green, IND);
    T = reshape(K4_grid(:,:,:,IND(i,:)), n(1),n(2),n(3), ndim, ndim);
    rhs = residual_func3D(T, Green);
    [SVect, ~] = pcg(lhs, rhs, 1e-10, 50);
    S(:,:,:,I(i), J(i)) = reshape(SVect, n_nodes, ndim, ndim);
end
K4_Macro = zeros(1, length(IND_K4));
for i = 1:length(IND_K4)
    idx = sub2ind([ndim, ndim], IND_K4(i,3), IND_K4(i,4));  SS = squeeze(S(:,:,:,idx));
    idx = sub2ind([ndim, ndim], IND_K4(i,1), IND_K4(i,2));
    K4_Macro(i) = h*(sum(K4(:,i)) + sum(K4(:,IND(idx,:)) .* SS(:,:), [1, 2])) / vol;
end

%% FIGURE: LAMINATE RVE
figure('Position', 0.6*get(0, 'DefaultFigurePosition'))
inv_transparent = 0.6;
plotcube([2/3, 2, 2], -1*ones(1,3), inv_transparent, [0.5, 0, 0]);
plotcube([2/3, 2, 2], [-1 + 2/3, -1, -1], inv_transparent, [0, 0.5, 0]);
plotcube([2/3, 2, 2], [-1 + 4/3, -1, -1], inv_transparent, [0, 0, 0.5]);
X = [-1, 1; -1 1; -1 1;  -1 1; -1 -1; -1 -1; -1 -1; -1 -1;  1 1; 1 1; 1 1; 1 1]';
Y = [-1, -1; 1 1; -1 -1; 1 1;   -1 -1; 1 1; -1 1; -1 1;  -1 -1; 1 1; -1 1; -1 1]';
Z = [-1, -1; 1 1; 1 1; -1 -1;   -1 1; -1 1; -1 -1; 1 1;  -1 1; -1 1; -1 -1; 1 1]';
line(X, Y, Z, 'LineWidth', 2, 'Color', 'k');
axis square

s = 1.2;
axis(reshape(repmat([-s, s], 3, 1)', 1, []))
ax = gca; ax.FontSize = 11; fs = 14;
xlabel('$X_1$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$X_2$', 'Interpreter', 'latex', 'FontSize', fs);
zlabel('$X_3$', 'Interpreter', 'latex', 'FontSize', fs);
if exist('./figures', 'dir') ~= 7, mkdir 'figures'; end
exportgraphics(gcf, './figures/RVE_3D_laminate.eps', 'ContentType', 'image');

%% ACCURACY PLOT FOR THE 3-PHASE LAMINATED MICROSTRUCTURE
% For the solution of deformation gradient
figure('Position', 0.6*get(0, 'DefaultFigurePosition'))
F_analytical = readmatrix('deformationGradient.csv');
F = reshape(F, [n, ndim, ndim]);
F1 = squeeze(F(fix(n(1)/6), 1, 1, :, 1));
F2 = squeeze(F(fix(n(1)/6)+n(1)/3, 1, 1, :, 1));
F3 = squeeze(F(fix(n(1)/6)+2*n(1)/3, 1, 1, :, 1));
F1_diff = F_analytical(:,1) - F1;
F2_diff = F_analytical(:,2) - F2;
F3_diff = F_analytical(:,3) - F3;
F_diff = vertcat(F1_diff(:), F2_diff(:), F3_diff(:));
plot(1:3, F1_diff, '-o', 1:3, F2_diff, '-o', 1:3, F3_diff, '-o', 'LineWidth', 1.2), grid on
xticks(1:3)
ax = gca; fs = 13;
ax.TickLabelInterpreter = 'latex';
xticklabels({'$p = 1$', '$p = 2$', '$p = 3$'})
xlabel('$F_{i1}^{(p)}$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$F_{i1}^{\mathrm{ref}} - F_{i1}^{\mathrm{FFT}}$', 'Interpreter', 'latex', 'FontSize', fs);
xlim([0.9, 3.1])
y_frame = max(F_diff) - min(F_diff);
ylim([min(F_diff) - 0.1 * y_frame, max(F_diff) + 0.1*y_frame])
legend('$F_{11}$', '$F_{21}$', '$F_{31}$', 'Interpreter', 'latex', ...
    'FontSize', 10, 'Location', 'north', 'box', 'off')
if exist('./figures', 'dir') ~= 7, mkdir 'figures'; end
exportgraphics(gcf, './figures/three_phase_solution_accuracy.eps', 'ContentType', 'image');

% For the solution of elastic moduli
figure('Position', 0.6*get(0, 'DefaultFigurePosition'))
K4 = readmatrix('ElasticModuli.csv');
K4_numerical = reshape(K4_Macro, [], 1);
K4_analytical = reshape(K4, [], 1);
relative_error = (K4_numerical - K4_analytical) ./ K4_analytical * 100;
plot( relative_error, 'k-*', 'LineWidth', 0.7); grid on
xticks(0:5:45);
xlabel('$K_{\mathrm{IND}}$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\big(K_\mathrm{IND}^\mathrm{num} - K_\mathrm{IND}^\mathrm{ref}\big)/K_\mathrm{IND}^\mathrm{ref} (\%)$', ...
    'Interpreter', 'latex', 'FontSize', fs)
xlim([0, 46])
if exist('./figures', 'dir') ~= 7, mkdir 'figures'; end
exportgraphics(gcf, './figures/tangent_moduli_comparison.eps', 'ContentType', 'image');

%% AUXILIARY FUNCTIONS
function v = residual_func3D(P, G)  % Compute right-hand side for the CG solver
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
function v = stiff_func3D(dF_vector, K4, G, IND) % Define linear function AFUN for the CG solver
n = zeros(1,3); [n(1), n(2), n(3), ~] = size(G); ndim = 3;
dF = reshape(dF_vector, [n, ndim, ndim]);    % unroll vector into tensor-matrix
dK = zeros([n, ndim, ndim]);             % tangent increment
for i = 1:length(IND), dK(:,:,:,i) = fftn( sum(K4(:,:,:, IND(i,:)) .* dF(:,:,:, :), 4)); end
v = zeros([n, ndim, ndim]);
for i = 1:ndim
    v(:,:,:,i,1) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[1, 2, 3]), 4) );
    v(:,:,:,i,2) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[2, 4, 5]), 4) );
    v(:,:,:,i,3) = ifftn( sum(squeeze(dK(:,:,:,i,:)).*G(:,:,:,[3, 5, 6]), 4) );
end
v = v(:);      % roll tensor into column vector.
end