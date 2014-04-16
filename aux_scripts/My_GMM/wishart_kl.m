function [D] = wishart_kl (B_q,B_p,alpha_q,alpha_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [D] = wishart_kl (B_q,B_p,alpha_q,alpha_p)
%
%   computes the divergence 
%                /
%      D(q||p) = | q(x)*log(q(x)/p(x)) dx
%               /
%   between two k-dimensional Wishart propability density for C  given
%   scale matrix B and shape parameter alpha
%
%            alpha
%         |B|              alpha-(k+1)/2
%   p(x)= ------------- |C|              exp (-tr(BC))
%         Gamma (alpha)
%              k
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4,
  error('Incorrect number of input arguements');
end;

if size(B_q)~=size(B_p),
  error('Distributions must have equal dimensions');
end;

K=size(B_p,1);

DBq=det(B_q);
DBp=det(B_p);

D=alpha_p*log(DBq/DBp)+alpha_q*(trace(B_p*inv(B_q))-K);
for k=1:K,
  D=D+gammaln(alpha_p+0.5-0.5*k)-gammaln(alpha_q+0.5-0.5*k);
  D=D+(alpha_q-alpha_p)*digamma(alpha_q+0.5-0.5*k);
end;

return;

