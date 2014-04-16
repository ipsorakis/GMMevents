function B = build_B_from_Y(DATA,Y)

individuals = unique(DATA(:,2));
N_individuals = length(individuals);
%N_data_points = size(Y,1);
K = size(Y,1);


B = zeros(K,N_individuals);

for i=1:N_individuals
    i_indices = DATA(:,2) == individuals(i);
    
    if sum(i_indices)~=1
        B(:,individuals(i)) = sum(Y(:,i_indices),2);
    else
        B(:,individuals(i)) = Y(:,i_indices);
    end
end

end