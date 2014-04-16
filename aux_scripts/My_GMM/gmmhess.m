%[H,L] = gmmhess(x,m,C,p,full_or_diag)
%
%x	data
%m	centres
%C	cov elements
%full_or_diag	'f' or 'd'
%

function [H,L] = gmmhess(x,m,C,p,full_or_diag)

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

H = K*log(N) - sum(log(Pj)) + 2*d*sum(log(sqrt(2)*Nj)) - 2*sum(Sj);

%a=2;
%s2_pop = mean(diag(cov(x)));
%Vo_p = 1;
%meano_p = ones(size(p,1),1)/K;
%Vo_m = a*a*s2_pop;
%Vo_s = s2_pop;

%ln_Pcov = K*d*log(2*s2_pop) - log(prod(1:K-1));
%ln_Pcov = 2*K*d*log(a*s2_pop) + K*log(Vo_p)
%Dtheta = sum(sum(D))/Vo_s + sum(sum(m.^2))/Vo_m + (p-meano_p)'*(p-meano_p)/Vo_p
%H = H + ln_Pcov;

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
