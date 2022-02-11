function [P, varargout] = neo_hookean(F, mu, beta, IND_K4)
[ndof, ndim, ~] = size(F);
J = F(:,1,1).*(F(:,2,2).*F(:,3,3) - F(:,2,3).*F(:,3,2)) ...
    + F(:,1,2).*(F(:,2,3).*F(:,3,1) - F(:,2,1).*F(:,3,3)) ...
    + F(:,1,3).*(F(:,2,1).*F(:,3,2) - F(:,2,2).*F(:,3,1));
JBeta = J.^(-beta);
invF = zeros(ndof, ndim, ndim);
invF(:,1,1) = F(:,2,2).*F(:,3,3) - F(:,2,3).*F(:,3,2); 
invF(:,1,2) = F(:,1,3).*F(:,3,2) - F(:,1,2).*F(:,3,3);
invF(:,1,3) = F(:,1,2).*F(:,2,3) - F(:,1,3).*F(:,2,2);
invF(:,2,1) = F(:,2,3).*F(:,3,1) - F(:,2,1).*F(:,3,3);
invF(:,2,2) = F(:,1,1).*F(:,3,3) - F(:,1,3).*F(:,3,1);
invF(:,2,3) = F(:,1,3).*F(:,2,1) - F(:,1,1).*F(:,2,3);
invF(:,3,1) = F(:,2,1).*F(:,3,2) - F(:,2,2).*F(:,3,1);
invF(:,3,2) = F(:,1,2).*F(:,3,1) - F(:,1,1).*F(:,3,2);
invF(:,3,3) = F(:,1,1).*F(:,2,2) - F(:,1,2).*F(:,2,1);
invF = invF./J;
P = zeros(size(F));    % 1st Piola-Kirchhoff stress.
P(:,1:9) = mu.*F(:,1:9) - mu.*invF(:,[1, 4, 7, 2, 5, 8, 3, 6, 9]).*JBeta;
if nargout > 1
    invF_x_invF = zeros(ndof, length(IND_K4));
    for q = 1:length(IND_K4)
        invF_x_invF(:,q) = invF(:,IND_K4(q,2), IND_K4(q,3)) .* invF(:,IND_K4(q,4),IND_K4(q,1));
    end
    K = zeros(ndof, ndim*ndim*(ndim*ndim+1)/2);
    for q = 1:length(IND_K4)
        [idx_bool, idx_loc] = ismember(IND_K4(q, [3, 2, 1, 4]), IND_K4, 'rows');
        if idx_bool == false
            [~, idx_loc] = ismember(IND_K4(q, [1, 4, 3, 2]), IND_K4, 'rows');
        end
        K(:,q) = mu.*JBeta.*(beta.*invF_x_invF(:,idx_loc) + invF_x_invF(:,q));
        if (IND_K4(q,1) == IND_K4(q,3)) && (IND_K4(q,2) == IND_K4(q,4))
            K(:,q) = K(:,q) + mu;
        end
    end
    varargout{1} = K;
end