function B = get_hard_incidence_matrix(P)
K = size(P,1);
N = size(P,2);

B = zeros(K,N);

for i=1:N
   k = find(P(:,i) == max(P(:,i)),1);
   B(k,i) = 1;
end
end