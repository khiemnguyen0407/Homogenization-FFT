function [P, varargout] = saint_venant(F, lambda, mu, IND_K4)
[ndof, ndim, ~] = size(F);
C = zeros(ndof, ndim, ndim);            % right Cauchy tensor
for i = 1:ndim
    for j = 1:i-1                       % compute only the symmetric part
        C(:,i,j) = sum(F(:,:,i) .* F(:,:,j), 2);
        C(:,j,i) = C(:,i,j);
    end
    C(:,i,i) = sum(F(:,:,i) .* F(:,:,i), 2);	% compute the diagonal components  
end     % Only the lower part of C is computed after this loop
E = zeros(ndof, ndim, ndim);            % Green-Lagrange tensor
E(:, [1, 5, 9]) = 0.5 * (C(:, [1, 5, 9]) - 1);
E(:, [2, 3, 6]) = 0.5 * C(:, [2, 3, 6]);
E(:, [4, 7, 8]) = E(:, [2, 3, 6]);      % for symmetric part of E

tr_green = sum(E(:, [1, 5, 9]), 2);     % trace of Green-Lagrange tensor
S = zeros(ndof, ndim, ndim);            % 2nd Pioloa-Kirchhoff stress
S(:,[1, 5, 9]) = lambda .* tr_green + 2*mu .* E(:,[1, 5, 9]);
S(:,[2, 3, 6]) = 2*mu .* E(:,[2, 3, 6]);
S(:,[4, 7, 8]) = S(:,[2, 3, 6]);        % for symmetric part of S
P = zeros(ndof, ndim, ndim);            % 1st Piola-Kirchhoff stress
for i = 1:ndim
    for j = 1:ndim
        P(:,i,j) = F(:,i,1).*S(:,1,j) + F(:,i,2).*S(:,2,j) + F(:,i,3).*S(:,3,j);
    end
end
if nargout > 0
    lambda_plus_2mu = lambda + 2*mu;
    lambda_plus_mu  = lambda + mu;
    F_x_F = zeros(ndof, length(IND_K4));
    for q = 1:length(IND_K4)
        F_x_F(:,q) = F(:,IND_K4(q,1),IND_K4(q,2)).*F(:,IND_K4(q,3),IND_K4(q,4)); 
    end
    K = zeros(ndof, ndim*ndim*(ndim*ndim+1)/2);
    % Compute the components K_{11xx}
    K(:,1) = lambda.*F_x_F(:,1) + mu.*(2*F_x_F(:,1) + F_x_F(:,10) + F_x_F(:,18)) + S(:,1,1);
    K(:,2) = lambda_plus_mu.*F_x_F(:,2) + S(:,2,1);
    K(:,3) = lambda_plus_mu.*F_x_F(:,3) + S(:,3,1);
    K(:,4) = lambda .* F_x_F(:,4) + mu.*(2*F_x_F(:,4) + F_x_F(:,13) + F_x_F(:,21));
    K(:,5) = mu .* F_x_F(:,12) + lambda .* F_x_F(:,5);
    K(:,6) = mu .* F_x_F(:,19) + lambda .* F_x_F(:,6);
    K(:,7) = lambda .* F_x_F(:,7) + mu.*(2*F_x_F(:,7) + F_x_F(:,16) + F_x_F(:,24));
    K(:,8) = mu .* F_x_F(:,15) + lambda .* F_x_F(:,8);
    K(:,9) = mu .* F_x_F(:,22) + lambda .* F_x_F(:,9);
    % Compute the components K_{12xx}
    K(:,10) = lambda.*F_x_F(:,10) + mu.*(F_x_F(:,1) + 2*F_x_F(:,10) + F_x_F(:,18)) + S(:,2,2);
    K(:,11) = lambda_plus_mu .* F_x_F(:,11) + S(:,3,2);
    K(:,12) = lambda.*F_x_F(:,12) + mu.*F_x_F(:,5);
    K(:,13) = lambda.*F_x_F(:,13) + mu.*(F_x_F(:,4) + 2*F_x_F(:,13) + F_x_F(:,21));
    K(:,14) = mu.*F_x_F(:,20) + lambda.*F_x_F(:,14);
    K(:,15) = lambda.*F_x_F(:,15) + mu.*F_x_F(:,8);
    K(:,16) = lambda.*F_x_F(:,16) + mu.*(F_x_F(:,7) + 2*F_x_F(:,16) + F_x_F(:,24));
    K(:,17) = mu.*F_x_F(:,23) + lambda.*F_x_F(:,17);
    % Compute the components K_{13xx}
    K(:,18) = lambda.*F_x_F(:,18) + mu.*(F_x_F(:,1) + F_x_F(:,10) + 2*F_x_F(:,18)) + S(:,3,3);
    K(:,19) = lambda.*F_x_F(:,19) + mu.*F_x_F(:,6);
    K(:,20) = lambda.*F_x_F(:,20) + mu.*F_x_F(:,14);
    K(:,21) = lambda_plus_2mu.*F_x_F(:,21) + mu.*(F_x_F(:,4) + F_x_F(:,13));
    K(:,22) = lambda.*F_x_F(:,22) + mu.*F_x_F(:,9);
    K(:,23) = lambda.*F_x_F(:,23) + mu.*F_x_F(:,17);
    K(:,24) = lambda_plus_2mu.*F_x_F(:,24) + mu.*(F_x_F(:,7) + F_x_F(:,16) );
    % Compute the components K_{21xx}
    K(:,25) = lambda.*F_x_F(:,25) + mu.*(2*F_x_F(:,25) + F_x_F(:,31) + F_x_F(:,36)) + S(:,1,1);
    K(:,26) = lambda_plus_mu.*F_x_F(:,26) + S(:,2,1);
    K(:,27) = lambda_plus_mu.*F_x_F(:,27) + S(:,3,1);
    K(:,28) = lambda.*F_x_F(:,28) + mu.*(2*F_x_F(:,28) + F_x_F(:,34) + F_x_F(:,39));
    K(:,29) = lambda.*F_x_F(:,29) + mu.*F_x_F(:,33);
    K(:,30) = lambda.*F_x_F(:,30) + mu.*F_x_F(:,37);
    % Compute the components K_{22xx}
    K(:,31) = lambda.*F_x_F(:,31) + mu.*(F_x_F(:,25) + 2*F_x_F(:,31) + F_x_F(:,36)) + S(:,2,2);
    K(:,32) = lambda_plus_mu .* F_x_F(:,32) + S(:,3,2);
    K(:,33) = lambda.*F_x_F(:,33) + mu.*F_x_F(:,29);
    K(:,34) = lambda.*F_x_F(:,34) + mu.*(F_x_F(:,28) + 2*F_x_F(:,34) + F_x_F(:,39));
    K(:,35) = lambda.*F_x_F(:,35) + mu.*F_x_F(:,38);
    % Compute the components K_{23xx}
    K(:,36) = lambda.*F_x_F(:,36) + mu.*(F_x_F(:,25) + F_x_F(:,31) + 2*F_x_F(:,36)) + S(:,3,3);
    K(:,37) = lambda.*F_x_F(:,37) + mu.*F_x_F(:,30);
    K(:,38) = lambda.*F_x_F(:,38) + mu.*F_x_F(:,35);
    K(:,39) = lambda_plus_2mu .* F_x_F(:,39) + mu.*(F_x_F(:,28) + F_x_F(:,34));
    % Compute the components K_{31xx}
    K(:,40) = lambda.*F_x_F(:,40) + mu.*(F_x_F(:,40) + F_x_F(:,43) + F_x_F(:,45)) + S(:,1,1);
    K(:,41) = lambda_plus_mu .* F_x_F(:,41) + S(:,2,1);
    K(:,42) = lambda_plus_mu .* F_x_F(:,42) + S(:,3,1);
    % Compute the components K_{32xx}
    K(:,43) = lambda.*F_x_F(:,43) + mu.*(F_x_F(:,40) + 2*F_x_F(:,43) + F_x_F(:,45)) + S(:,2,2);
    K(:,44) = lambda_plus_mu .* F_x_F(:,44) + S(:,3,2);
    % Compute the components K_{33xx}
    K(:,45) = lambda.*F_x_F(:,45) + mu.*(F_x_F(:,40) + F_x_F(:,43) + 2*F_x_F(:,45)) + S(:,3,3);
    varargout{1} = K;
end