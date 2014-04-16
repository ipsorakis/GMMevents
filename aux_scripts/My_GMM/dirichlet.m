function [px] = dirichlet(x,alpha,logoption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [px] = dirichlet (x,alpha,logoption)
%
%   computes k-dimensional Dirichlet propability density for x
%   given prior counts alpha, with alpha_i >0 for all i=1:k
%
%            
%         Gamma(sum_{i=1}^k alpha_i)     alpha_1-1        alpha_k-1 
%   p(x)= --------------------------- x_1         .... x_k
%         prod_{i=1}^k Gamma(alpha_i)
%              
%
%   if logoption is set to one, the log-propability is returned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,
  logoption=0;
end;
epsscale=10;

alpha=alpha(:)';
k=length(alpha);
if (size(x,2)~=k),
  x=x';
end;
[N,k2]=size(x);
if k2~=k, error('x and alpha must be equally long'); end;
if any(alpha<=0) error('alphas must be > 0 '); end;
if any(abs(1-sum(x,2))>epsscale*eps) 
  error('Distributions over x must be complete (sum(x)=1)'); 
end;

normconst=gammaln(sum(alpha))-sum(gammaln(alpha));

if any(x(:)==0)
  px=ones(N,1);
  for i=1:k
    px(:)=px.*x(:,i).^(alpha(i)-1);
  end;
  if logoption
    px=normconst+log(px);
  else
    px=exp(normconst)*px;
  end
else
  px=zeros(N,1);
  for i=1:k
    px(:)=px+(alpha(i)-1)*log(x(:,i));
  end;
  if logoption
    px=normconst+px;
  else 
    px=exp(normconst+px);
  end;
end;


