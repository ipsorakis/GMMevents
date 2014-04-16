function emc = em_circ(a,K,w);

%emc = em_circ(a,K,w);
%em_circ returns K wrapped normals on circular data a
%optional weighting w>0 (same size as a)

tol = 0.01; % nominal tolerance
MaxIter = 100; % nominal maximum iterations
UniWeight = 0.1; % relative prior weight of Uniform disn

if (sum(a < 0) > 0) % there is -ve data
    disp('correcting to range 0,2*pi');
    a = a+pi;
end;

if (max(a)>2*pi) % something's wrong
    error('data is out of corrected range 0,2*pi');
end;

if (nargin < 3) %no weights included
    w = ones(size(a));
end;

% initialise
pxgu = 1/(2*pi); % uniform in interval 0,2*pi
km = km_circ(a,K);
m = km.m;
C = km.C;
pk = km.pk*(1-UniWeight); % assume UniWeight proportion of data from Uniform
%m = rand(1,K)*2*pi; % random in interval 0,2*pi
%m = linspace(2*pi/K,2*pi,K); % even in interval 0,2*pi
%C = 4*ones(K,1)*(2*pi/K)^2; % std as twice intermean distance to start with
%pk = 0.8*ones(K,1)/(K+1); % flat priors to start with, with single uniform disn
pu = 1-sum(pk);
pp = zeros(length(a),K); % posteriors

Qnew = 0; % to start with

for iter = 1:MaxIter,
    if (rem(iter,10)==0), fprintf('.');end; % PRD
    Qold = Qnew;
    % E-step first to get posteriors
    for k=1:K
        z = a;
        e1 = m(k) - z;
        e2 = m(k) - (z+2*pi);
        d(:,k) = min(e1.^2, e2.^2);
        f = find(e2.^2 < e1.^2); % all those that need a wrapping
        z(f) = z(f) + 2*pi; % wrap those that need it
        pxgk = exp(-0.5*d(:,k)/C(k))/sqrt(2*pi*C(k));
        pp(:,k) = pxgk*pk(k);
    end;
    pugx = pxgu*pu;
    sum_pp = sum(pp,2)+pugx;
    Qnew = sum(log(sum_pp));
    pugx = pugx/sum_pp;
    
    % M-step to update m, C & pk & uniform disn
    for k=1:K
        pp(:,k) = pp(:,k)./sum_pp; % normalis the pps
        z = a;
        e1 = m(k) - z;
        e2 = m(k) - (z+2*pi);
        d(:,k) = min(e1.^2, e2.^2);
        f = find(e2.^2 < e1.^2); % all those that need a wrapping
        z(f) = z(f) + 2*pi; % wrap those that need it
        m(k) = sum(w.*pp(:,k).*z)/sum(w.*pp(:,k));
        m(k) = mod(m(k),2*pi);
        C(k) = sum( ((z-m(k)).^2).*pp(:,k).*w)/sum(w.*pp(:,k));
        pk(k) = mean(pp(:,k));
    end;
    pu = 1-sum(pk);
    if ( abs(Qnew-Qold)/abs(Qnew) < tol)
        StrOut=sprintf('Tolerance level reached after %d iterations\n',iter);
        disp(StrOut);
        break;
    end;
end; % iterations
fprintf('\n');

[dum,emc.k_star] = max([pp,pugx']');
emc.m = m;
emc.C = C;
emc.pk = pk;
emc.pu = pu;
emc.pp = pp;
[dum,ev1] = getposts(m,C,pk,a,'d');
[dum,ev2] = getposts(m,C,pk,a+2*pi,'d');
emc.ev = max(ev1,ev2);

