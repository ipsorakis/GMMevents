function emc = em_circ_multid(a,K, m, C, pk);
%em_circ returns K wrapped normals on circular data a

tol = 0.001; % nominal tolerance
MaxIter = 100; % nominal maximum iterations
UniWeight = 0.1; % relative prior weight of Uniform disn

[I, J]= find(a <0);
for i=1:length(I),
  a(I(i), J(i)) = a(I(i), J(i))+ 2*pi;
end

[I, J]= find(a > 2*pi);
for i=1:length(I),
  a(I(i), J(i)) = a(I(i), J(i)) - 2*pi;
end

if (sum(a < 0) > 0) % there is -ve data
    disp('correcting to range 0,2*pi');
    a = a+pi;
end;

if (max(a)>2*pi) % something's wrong
    error('data is out of corrected range 0,2*pi');
end;

% initialise
[L, D] = size(a);
pxgu = 1/(2*pi);
pk = pk*(1-UniWeight); % assume UniWeight proportion of data from Uniform
pu = 1-sum(pk);

pp = zeros(length(a),K); % posteriors
Qnew = 0; % to start with

for k=1:K
    sigma = sbmatout(C, D, k);
    if (det(sigma)<0.00001)
        b = 0.0005;  
        sigma  = sigma + b*eye(D);
        C = sbmatin(C, sigma, k);
    end
end




disp('em_circ');
for iter = 1:MaxIter,
    %if (rem(iter,10)==0), fprintf('.');end; % PRD
    fprintf('.');
    Qold = Qnew;
 %   disp('here');
    % E-step first to get posteriors
    disp('---------------');
    for k=1:K
        z = a;
        z = data_circ(z, m(k,:)); % wrap the data 
        sigma = sbmatout(C, D, k); %get the covariance matrix for the cluster
        mm = repmat(m(k,:), L,1);
        d = z-mm;
         disp(['det sig=', num2str(det(sigma)) ])
         disp(['pk=', num2str(pk(k))])

         for i=1:L
         pxgk(i) = exp(-0.5*d(i,:)*(inv(sigma))*d(i,:)')/ ((2*pi)^(D/2)*sqrt(det(sigma)));% make this the multivariate equation exp(-0.5*d(:,k)/C(k))/sqrt(2*pi*C(k));
         end
         if (pxgk>1.0)
            pxgk 
         end
        pp(:,k) = pxgk*pk(k);
    end;
    pugx = pxgu*pu;
    sum_pp = sum(pp,2)+pugx;
    KVL = find(sum_pp > 1.0);
    length(KVL)
    Qnew = sum(log(sum_pp));
    pugx = pugx/sum_pp;
    
    % M-step to update m, C & pk & uniform disn
    for k=1:K
        pp(:,k) = pp(:,k)./sum_pp; % normalis the pps
        z = a;
        z = data_circ(z, m(k,:)); % wrap the data 
        pprep = repmat(pp(:,k), 1, D);
        m(k,:) = sum(pprep.*z)/sum(pp(:,k)); % will this give a vector
        m(k,:) = mod(m(k,:),2*pi);
        %C(k) = sum( ((z-m(k)).^2).*pp(:,k))/sum(pp(:,k));
        mm = repmat(m(k,:), L,1);
        d = z-mm;
        sigma = zeros(D,D);
        for i=1:L,
            sigma = sigma + (d(i,:)'*d(i,:)*pp(i,k));
        end
        sigma = sigma/sum(pp(:,k));
        
        if (det(sigma)<0.00001)
          b = 0.0005;  
          sigma  = sigma + b*eye(D);
        end
            
        %disp('size sigma:');
        %size(sigma)
        C = sbmatin(C, sigma, k);
        pk(k) = mean(pp(:,k));
    end;
    %pu = 1-sum(pk);
    if ( abs(Qnew-Qold)/abs(Qnew) < tol)
        StrOut=sprintf('Tolerance level reached after %d iterations\n',iter);
        disp(StrOut);
        break;
    end;
end; % iterations
fprintf('\n');

[dum,emc.k_star] = max([pp,pugx']');
%[dum,emc.k_star] = max(pp');
emc.m = m;
emc.C = C;
emc.pk = pk;
emc.pu = pu;
emc.pp = pp;
emc.logLikelihood = Qnew;
%[dum,ev1] = getposts(m',C',pk,a,'d');
%[dum,ev2] = getposts(m',C',pk,a+2*pi,'d');
%[dum,ev3] = getposts(m',C',pk,a-2*pi,'d');
%emc.ev = max([ev1,ev2, ev3]);

