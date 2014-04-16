function [x]=wishartrnd(alpha,C,N)
%
%  x=WISHARTRND(alpha,C,N)
%
%  samples N-times from an K-dimensional Wishart distribution with kxk
%  scale matrix C and  alpha degrees for freedom. Dimensionality is implied
%  in the scale matrix
%
%  e.g: C=[1 .7;0.7 1];
%       alpha=10;
%       x=wishartrnd(alpha,C,300);

% We use the method from p. 480 of "Bayesian Data Analysis", Gelman et al.

alpha=ceil(alpha);
ndim=size(C,1);

if ndim>alpha,
  error('D.o.f. must be larger than dimensionality');
end;
if det(C)<=0
  error('Scale matrix must be positive definite');
end;

m=zeros(ndim,1);
y=sampgauss(m,C,alpha);
x=y*y';

for i=2:N,
  y=sampgauss(zeros(ndim,1),C,alpha);
  x(:,:,i)=y*y';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=sampgauss(m,C,N)
%
%  x=SAMPGAUSS(m,C,N)
%
%  samples N-times from an multi-dimensional gaussian distribution 
%  with covariance matrix C and mean m. Dimensionality is implied
%  in the mean vector
%
%  e.g: C=[1 .7;0.7 1];
%       m=[0;0];
%       x=sampgauss(m,C,300);
m=m(:);
ndim=size(C,1);
if size(C,2)~=ndim,
  x=[];
  error('Wrong specification calling sampgauss');
end

if ndim==1,
   x=m+C*randn(1,N);
   return;
end;

% check determinant of covariance matrix
%if det(C)>1 | det(C)<=0, error('Covariance matrix determinant must be
%0< |C| <= 1'); end;

% generate zero mean/unit variance samples
e=randn(ndim,N);
% make sure they are unit variance;
s=std(e');
for i=1:ndim,
   e(i,:)=e(i,:)./s(i);
end;

% decompose cov-matrix
S=chol(inv(C));
x=inv(S)*e+m(:,ones(N,1));
