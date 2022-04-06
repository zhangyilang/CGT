function [x_mols, x_supp, iter_count] = MOLS_cK(y, Phi, c, K, rt, S)
%MOLS_CK_PRUNE

[nRows]     	= size(Phi,1);
residual_prev	= y;
supp			= [];
iter_count		= 0; 
zero_threshold  = rt; 
max_itr			= min(c * K, floor(nRows/S)); % maximum iteration number

while (norm(residual_prev) > zero_threshold && iter_count < max_itr)
    iter_count	= iter_count + 1;

    if 0 == size(supp) 
        Phi1    = Phi;    
    else
        Phi1    = (eye(size(Phi,1)) - Phi(:,supp)*pinv(Phi(:,supp))) * Phi;
        Phi1    = Phi1*diag(1./sparse(sqrt(sum(Phi1.^2,1)))); 
    end

    corr        = abs(Phi1'*residual_prev);  
    corr(supp)  = 0;     
    [~, supp_idx]	    = sort(corr, 'descend');
    supp_n				= union(supp, supp_idx(1:S));

    if isempty( intersect(supp,supp_idx(1:S)) ) && (length(supp_n) <= nRows )
        x_hat			= iplsp_EstLS(y, Phi(:,supp_n));
        residual_prev	= y - Phi(:,supp_n)*x_hat;
        supp	        = supp_n;
    else
        break; % in this case, supp os not overwritten
    end
end

x_supp			         = supp;
x_mols(supp)	         = iplsp_EstLS(y, Phi(:,x_supp));

% Re-estimate the value in x_mols via an LS projection with respect to K largest magnitudes in x_hat 
% [~, supp_idx]            = sort(abs(x_mols), 'descend');
% x_supp                   = supp_idx(1:min(c * K, length(supp_idx)));
% x_mols                   = zeros(size(Phi,2), 1);
% x_mols(x_supp)           = iplsp_EstLS(y, Phi(:,x_supp));

end

