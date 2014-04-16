function B = build_pulse_bipartite_graph(DATA,Y)

individuals = unique(DATA(:,2));
N_individuals = length(individuals);
%N_data_points = size(Y,1);
K = size(Y,2);


B = zeros(N_individuals,K);

for i=1:N_individuals
    i_indices = DATA(:,2) == i;
    
    n_i = sum(i_indices);
    
    if n_i > 1
        B(i,:) = (1/n_i) * sum(Y(i_indices,:));
    else
        B(i,:) = Y(i_indices,:);
    end
end

end