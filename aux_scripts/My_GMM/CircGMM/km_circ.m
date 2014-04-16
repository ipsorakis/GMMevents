function kmc = km_circ_multi(a,K);

%kmc = km_circ_multi(a,K);
%km_circ_multi returns K means on multi-d circular data a with K means

tol = 0; % nominal tolerance
MaxIter = 100; % nominal maximum iterations

[L,D] = size(a);
if (D>L)
     a = a'; %swap round
end;

if (sum(sum(a < 0)) > 0) % there is -ve data
    disp('correcting to range 0,2*pi');
    a = a+pi;
end;

if (max(max(a))>2*pi) % something's wrong
    error('data is out of corrected range 0,2*pi');
end;

% now can start the K-means

%m = rand(1,K)*2*pi; % random in interval 0,2*pi
m = linspace(2*pi/K,2*pi,K)'; % even in interval 0,2*pi
m = repmat(m,1,D);
C=[];

for iter = 1:MaxIter,
    old_m = m;
    if (rem(iter,10)==0), fprintf('.');end; % PRD

    [se,k_star] = assign_km_labels(a,m);

    for k=1:K
        f = find(k_star==k); % all the relevant vectors
        z = a(f,:);
        mm = repmat(m(k,:),length(z),1);
        z = data_circ(z,m(k,:));
        m(k,:) = mean(z);
        C = sbmatin(C,cov(z-mm),k);
    end;
    
    if ( norm(m-old_m,2)/norm(m,2) < tol)
        StrOut=sprintf('Tolerance level reached after %d iterations\n',iter);
        disp(StrOut);
        break;
    end;
end; % iterations
fprintf('\n');

[se,k_star] = assign_km_labels(a,m);
m = mod(m,2*pi);

for k=1:K
   p(k) = sum(k_star==k);
end;
p = p/length(k_star);

kmc.pk = p;
kmc.m = m;
kmc.C = C;
kmc.k_star = k_star;
kmc.energy = sum(se);
[dum,ev1] = getposts(m,C,p,a,'d');
[dum,ev2] = getposts(m,C,p,a+2*pi,'d');
kmc.ev = max(ev1,ev2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [se,k_star] = assign_km_labels(a,m);

[K,D] = size(m);
[L,D] = size(a);

for k=1:K
  mm = repmat(m(k,:),L,1);
  a = data_circ(a,m(k,:));
  d(:,k) = sum((a-mm).^2,2);
end;
[se,k_star] = min(d');
% end of function assign_km

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
