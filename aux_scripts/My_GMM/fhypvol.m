%[fhv] = fhypvol(C,K,full_or_diag);
%
%K	number of kernels
%C	cov elements
%full_or diag	'f' or 'd'
%


function [fhv] = fhypvol(C,K,full_or_diag);

if full_or_diag == 'd'
  D = C;
else
  D = diagcov(C,K);
end;

fhv = sum( sqrt(prod(D')) );

return;
