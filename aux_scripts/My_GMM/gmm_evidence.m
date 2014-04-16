%[E] = gmm_evidence(th,x,K)
%
%x	data
%th	model paremeters
%K	number of kernels
%E	evidence term

function [E] = gmm_evidence(th,x,K)

d = size(x,2);
N = size(x,1);

m = reshape(th(1 : K*d),K,d);
C = reshape(th(K*d+1 : 2*K*d),K,d);
%Pj = reshape(th(2*K*d+1 : 2*K*d+K),K,1);
Pj = ones(K,1)/K;

Nj = Pj.*N;
D = C;
Sj = sum(log(D)')';

H = K*log(N) - sum(log(Pj)) + 2*d*sum(log(sqrt(2)*Nj)) - 2*sum(Sj);

%a=2;
%s2_pop = mean(diag(cov(x)));
%Vo_p = 1;
%meano_p = ones(size(p,1),1)/K;
%Vo_m = a*a*s2_pop;
%Vo_s = s2_pop;

%ln_Pcov = K*d*log(2*s2_pop) - log(prod(1:K-1));
%ln_Pcov = 2*K*d*log(a*s2_pop) + K*log(Vo_p);
%Dtheta = sum(sum(D))/Vo_s + sum(sum(m.^2))/Vo_m + (p-meano_p)'*(p-meano_p)/Vo_p;
%H = H + ln_Pcov + Dtheta;

for j=1:K
  F = (ones(d,1)*C(j,:)).*eye(d);
  Pjx(j,:) = gaussres(x,m(j,:),F);
end;

L = sum( log(Pj'*Pjx) );
E = -L + 0.5*H;

return;
