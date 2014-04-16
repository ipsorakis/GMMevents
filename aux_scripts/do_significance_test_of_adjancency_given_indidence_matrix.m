function [W_significant W_null_mean W_null_std P_values P_significant Wall] = do_significance_test_of_adjancency_given_indidence_matrix(W,B,a_threshold,number_of_randomisations)

[K N] = size(B);
%P_values = zeros(N);

P_values_left = zeros(N);
P_values_right = zeros(N);


Wall = zeros(N^2,number_of_randomisations);
%Wall = zeros(N,N,number_of_randomisations);

parfor rand_index = 1:number_of_randomisations
    %% randomise B
    %    Baux = B;
    %    Bnull = Baux(:);
    %    Bnull = Bnull(randperm(N),:);
    Bnull = B;
    for n=1:N
        visitation_profile = Bnull(:,n);
        occurences_n = sum(visitation_profile);
        pi_k = visitation_profile / occurences_n;
        pi_k = pi_k(randperm(K));
        
        Bnull(:,n) = mnrnd(occurences_n,pi_k)';
    end    
    
    %% find Wnull
    Wnull = get_coocurences_in_bipartite_graph(Bnull);
    
    Wall(:,rand_index) = Wnull(:);
end

for k=1:N^2
    j = floor((k-1)/N) + 1;
    i = k - (j-1)*N;
    
    h = histc(Wall(k,:),0:max([Wall(k,:) W(i,j)]));
    h = h/sum(h);
    
    P_values_left(i,j) = sum(h(1:W(i,j)+1));
    P_values_right(i,j) = sum(h(W(i,j)+1:end));
    
    %P_values(i,j) = do_simple_p_test(h,W(i,j),1);
end

P_significant = zeros(N);
P_values = nan(N);
for i=1:N-1
    for j=i:N
        %P_significant(i,j) = (P_values_left(i,j)<a_threshold) || (P_values_right(i,j)<a_threshold);
        P_significant(i,j) = (P_values_right(i,j)<a_threshold);
        %P_values(i,j) = min(P_values_left(i,j),P_values_right(i,j));
        P_values(i,j) = P_values_right(i,j);
    end
end
P_significant = P_significant + P_significant';

W_significant = W.*P_significant;

A = W~=0;
P_significant = P_significant.*A;

%original_number_of_links = sum(get_triu_vector(A));
%significant_number_of_links = sum(get_triu_vector(P_significant));
%prune_factor = 1 - significant_number_of_links/original_number_of_links;

W_null_mean = mean(Wall,2);
W_null_mean = reshape(W_null_mean,N,N);

W_null_std = std(Wall');
W_null_std = reshape(W_null_std,N,N);
end