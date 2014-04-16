function T = get_T_matrix(Y,DATA)

K = size(Y,2);

T = zeros(K,2);

for k=1:K
    first_obs = find(Y(:,k)~=0,1,'first');
    last_obs = find(Y(:,k)~=0,1,'last');
    
    T(k,1) = DATA(first_obs,1);
    T(k,2) = DATA(last_obs,1);
end

end