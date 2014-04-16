function [px] = wishart(C,B,alpha,logoption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [px] = wishart (C,B,alpha,logoption)
%
%   computes k-dimensional Wishart propability density for C  given
%   scale matrix B and shape parameter alpha
%
%            alpha
%         |B|              alpha-(k+1)/2
%   p(x)= ------------- |C|              exp (-tr(BC))
%         Gamma (alpha)
%              k
%
%
%   if logoption is set to one, the log-propability is returned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4,
  logoption=0;
end;
epsscale=100;

k=size(B,1);

% testing C for squareness, symmetry; testing positive-definiteness in
% loop 
if size(C,1)~=k | size(C,2)~=k, error('C must be square'); end;
N=size(C,3);				% length of C
Ct=permute(C,[2 1 3]);			% transpose C
if any(any(abs(reshape(C,k*k,N)-reshape(Ct,k*k,N))>epsscale*eps,1)==1)
  error('C must be symmetric');
end;

% testing B for squareness, singularity and symmetry
if size(B,1)~=k | size(B,2)~=k, error('B must be square'); end;
Bt=B';
if any(B(:)~=Bt(:)), error('B must be symmetric'); end;
if rank(B)~=k, error('B must be non-singular'); end;

% testing for alpha 2*alpha > k-1
if 2*alpha<=(k-1) error('alpha must be greater than (k-1)/2'); end;


gammak=(k*(k-1)/4)*log(pi)+ sum(gammaln(alpha+.5-0.5*(1:k)));
normconst=alpha*log(det(B))-gammak;

px=zeros(N,1);
for i=1:N,
  detC=det(C(:,:,i));
  if detC<=0 error('C must be positive-definite'); end;
  px(i) = (alpha-0.5*(k+1))*log(detC) - trace(B*C(:,:,i));
end;

if logoption
  px=px+normconst;
else
  px=exp(px+normconst);
end;

