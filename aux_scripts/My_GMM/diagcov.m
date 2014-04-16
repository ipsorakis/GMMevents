%[D] = diagcov(FC,N)
%
%FC : full covariance matrix array
%N  : number of covs in array
%D  : diagonal eigenvalue array

function [D] = diagcov(FC,N)

dim = size(FC,2);
D = zeros(N,dim);

for n = 1:N
	A = sbmatout(FC,dim,n);
	D(n,:) = eig(A)';
end;				% for

return;


