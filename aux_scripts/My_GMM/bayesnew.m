%[H,L] = bayesnew(x,m,C,p,full_or_diag)
%
%x	data
%m	centres
%C	cov elements
%full_or_diag	'f' or 'd'
%

function [H,L] = bayesnew(x,m,C,p,full_or_diag)

d = size(x,2);
N = size(x,1);
K = size(m,1); 

Pj = p;
Nj = Pj.*N;
if full_or_diag == 'd'
  D = C;
else
  D = diagcov(C,K);
end;
Sj = sum(log(D)')';

[pp,ev] = getposts(m,C,p,x,full_or_diag);
ppn = pp./(ones(length(pp),1)*Pj');

if (K > 1)
  PT = sum( (pp(:,1:K-1) - pp(:,K)*ones(1,K-1)).^2 );
  PT = sum(log(PT));
else
  PT = log(N);
end;

%% original term 
%% H = K*log(N) - sum(log(Pj)) + 2*d*sum(log(sqrt(2)*Nj)) - 2*sum(Sj);
H = PT  + 2*d*sum(log(sqrt(2)*Nj)) - 2*sum(Sj);
H = H*0.5;

% evaluate the number of free parameters in the model (Np)

%%if full_or_diag == 'd'
  Np = d*K + d*K + K;
%%else
%%  Np = d*K + 0.5*d*(d+1)*K + K;
%%end;

H = H - (Np/2)*log(2*pi);

% now have 1/2 log |H| - p/2 log(2 pi)
% put in prior terms ( log P(theta) )

a=1;
b=1;
s2_pop = mean(diag(cov(x)));
ln_P = -K*d*log(2*a*b*s2_pop) + log(prod(1:K-1));

H = H - ln_P;

for j=1:K
  if full_or_diag == 'd'
    F = (ones(d,1)*C(j,:)).*eye(d);
  else
    F = sbmatout(C,d,j);
  end;
  Pjx(j,:) = gaussres(x,m(j,:),F);
end;

L = sum( log(Pj'*Pjx) );

return;
