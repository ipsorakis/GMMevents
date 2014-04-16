%[mml] = gmmmml(x,m,C,p,full_or_diag)
%
%x	data
%m	centres
%C	cov elements
%full_or_diag	'f' or 'd'
%

function [mml] = gmmmml(x,m,C,p,full_or_diag)

Kappa = [.08333,.080188,.077875,.07609,.07465,.07347,.07248,.07163];

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

% Num params estimation : centres, cov, priors (in that order)

%%if full_or_diag == 'd'
  Np = d*K + d*K + K;
%%else
%%  Np = d*K + 0.5*d*(d+1)*K + K;
%%end;

Fish_term = 0;

for j=1:K
  if full_or_diag == 'd'
    F = (ones(d,1)*C(j,:)).*eye(d);
  else
    F = sbmatout(C,d,j);
  end;
  Pjx(j,:) = gaussres(x,m(j,:),F);
  Fish_term = Fish_term + log(sqrt(2)*Nj(j)/prod(D(j,:)));
end;

L = sum( log(Pj'*Pjx) );
s2_pop = mean(diag(cov(x)));

if (Np > 8)
  lattice = Kappa(8);
else
  lattice = Kappa(Np);
end;
% first the prior terms
mml = K*d*log(2*s2_pop) - log(prod(1:K-1)) + 0.5*Np*log(lattice) - log(prod(1:K));
% now the expected Fisher Inf. terms
mml = mml + Fish_term + 0.5*log(N) - 0.5*sum(log(Pj));
% now likelihood terms
mml = mml - L + 0.5*Np;

return;
