function [D] = dirichlet_kl(alpha_q,alpha_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [D] = dirichlet_kl (alpha_q,alpha_p)
%
%   computes the divergence 
%                /
%      D(q||p) = | q(x)*log(q(x)/p(x)) dx
%               /
%   between two k-dimensional Dirichlet propability densities, where the pdf
%   of a dirichlet is given by  
%            
%         Gamma(sum_{i=1}^k alpha_i)     alpha_1-1        alpha_k-1 
%   p(x)= --------------------------- x_1         .... x_k
%         prod_{i=1}^k Gamma(alpha_i)
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,
  error('Incorrect number of input arguements');
end;

if length(alpha_q)~=length(alpha_p),
  error('Distributions must have equal dimensions');
end;
K=length(alpha_q);

aqtot=sum(alpha_q);
aptot=sum(alpha_p);
Psiqtot=digamma(aqtot);

D=gammaln(aqtot)-gammaln(aptot);
for k=1:K,
  D=D+gammaln(alpha_p(k))-gammaln(alpha_q(k))+...
      +(alpha_q(k)-alpha_p(k))*(digamma(alpha_q(k))-Psiqtot);

end;


