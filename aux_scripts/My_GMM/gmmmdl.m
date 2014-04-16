%[mdl] = gmmmdl(x,m,C,p,full_or_diag)
%
%x	data
%m	centres
%C	cov elements
%full_or_diag	'f' or 'd'
%

function [mdl] = gmmmdl(x,m,C,p,full_or_diag)

d = size(x,2);
N = size(x,1);
K = size(m,1); 

Pj = p';
Nj = Pj.*N;
if full_or_diag == 'd'
  D = C;
else
  D = diagcov(C,K);
end;

% Num params estimation : centres, cov, priors (in that order)

%%if full_or_diag == 'd'
  Np = d*K + d*K + K;
%%else
%%  Np = d*K + 0.5*d*(d+1)*K + K;
%%end;

for j=1:K
  if full_or_diag == 'd'
    F = (ones(d,1)*C(j,:)).*eye(d);
  else
    F = sbmatout(C,d,j);
  end;
  Pjx(j,:) = gaussres(x,m(j,:),F);
end;

L = sum( log(Pj'*Pjx) );

mdl = -L + 0.5*Np*log(N);

return;
