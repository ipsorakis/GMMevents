% [m,C,p] = gmmopt(x,K,m,C,P)
%
% optimises using BFGS 
%
% Inputs
%	x	data set
%	K	Number of kernels
% Outputs
%	m,C,p	optimised parameters

function [m,C,p] = gmmopt(x,K,m,C,p)

% sets up constants etc
[N,d] = size(x);
options = foptions;
options(1) = 1;		% show full display
options(2)=  1e-3;      % Tolerance
options(3) = 0.0;	% nmse criterion
options(14) = 500;	% number of iterations for each pass
options(7) = 1;		% line search routine
%%options(9) = 1;	% check my gradients

%%th = [randn(K*d,1);ones(K*d,1);ones(K,1)/K];
th(1:K*d) = reshape(m,K*d,1);
th(K*d + 1: 2*K*d) = reshape(C,K*d,1);
th(2*K*d + 1 : 2*K*d + K) = p;
th = bfgs_tol('gmm_evidence',th,options,[],x,K);

m = reshape(th(1 : K*d),K,d);
C = reshape(th(K*d+1 : 2*K*d),K,d);
p = reshape(th(2*K*d+1 : 2*K*d+K),K,1);

return;
