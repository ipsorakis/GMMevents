function W = get_coocurences_in_bipartite_graph(B)
% assumes [K x N] incidence matrix - K clusters N individuals
N = size(B,2);
K = size(B,1);

W = zeros(N);

for i=1:N-1
    for j=i+1:N
        
        visitation_profile_i = B(:,i);
        visitation_profile_j = B(:,j);
        
        % to vectorise
        for k=1:K
            W(i,j) = W(i,j) + min(visitation_profile_i(k),visitation_profile_j(k));

        end
    end
end

W = W + W';
end