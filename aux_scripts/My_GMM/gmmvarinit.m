function [gmm] = gmmvarinit(data,K,options);
% initialises the gaussian mixture model for Variational GMM algorithm
  MIN_COV=eps;

  
  gmm.K=K;
  [N,ndim]=size(data);

  midscale=median(data)';
  drange=max(data)-min(data);	
  % educated guess with scaling
  
  % define P-priors
  defgmmpriors=struct('Dir_alpha',[],'Norm_Mu',[],'Norm_Cov', ...
		      [],'Norm_Prec',[],'Wish_B',[],'Wish_iB',[],...
		      'Wish_alpha',[],'Wish_k',[]);
  
  defgmmpriors.Dir_alpha=ones(1,K);
  defgmmpriors.Norm_Mu=midscale;
  defgmmpriors.Norm_Cov=diag(drange.^2);
  defgmmpriors.Norm_Prec=inv(defgmmpriors.Norm_Cov);
  defgmmpriors.Wish_B=diag(drange);
  defgmmpriors.Wish_iB=inv(defgmmpriors.Wish_B);
  defgmmpriors.Wish_alpha=ndim+1;
  defgmmpriors.Wish_k=ndim;
  
  % assigning default P-priors 
  if ~isfield(options,'priors'),
    gmm.priors=defgmmpriors;
  else
    % priors not specified are set to default
    gmmpriorlist=fieldnames(defgmmpriors);
    fldname=fieldnames(gmm.priors);
    misfldname=find(~ismember(gmmpriorlist,fldname));
    for i=1:length(misfldname),
      priorval=getfield(defgmmpriors,gmmpriorlist{i});
      gmm.priors=setfield(gmm.priors,gmmpriorlist{i},priorval);
    end;
  end;

  % initialise posteriors
  switch options.init
   case 'rand'
    % sample mean from data
    ndx=floor(rand(1,K)*N+1);
    Mu=data(ndx,:)';
    %Mu=gaussrnd(gmm.priors.Norm_Mu,gmm.priors.Norm_Cov,K);
    % sample precision
    alpha=gmm.priors.Wish_alpha*2;
    Prec=wishartrnd(alpha,gmm.priors.Wish_iB,K);
    % sample weights
    kappa=gmm.priors.Dir_alpha;
    % assign values
    for k=1:K,
      % mean posterior
      gmm.post(k).Norm_Prec=gmm.priors.Norm_Prec;
      gmm.post(k).Norm_Cov=gmm.priors.Norm_Cov;
      gmm.post(k).Norm_Mu=Mu(:,k);
      % covariance posterior
      gmm.post(k).Wish_alpha=alpha;
      gmm.post(k).Wish_iB=squeeze(Prec(:,:,k))/alpha; %gmm.priors.Wish_B ?
      gmm.post(k).Wish_B=inv(gmm.post(k).Wish_iB);
      % weights posterior
      gmm.post(k).Dir_alpha=kappa(k);
    end;
   case 'conditional'
    [centre,Sigma,weight,Nm]=kminit(data,K,options);
    for k=1:K,
      Cov=squeeze(Sigma(:,:,k));
      Prec=inv(Cov);
      % mean posterior
      gmm.post(k).Norm_Prec=Nm(k)*Prec+gmm.priors.Norm_Prec;
      gmm.post(k).Norm_Cov=inv(gmm.post(k).Norm_Prec);
      gmm.post(k).Norm_Mu=gmm.post(k).Norm_Cov *...
	  (Nm(k)*Prec*centre(:,k)+...
	   gmm.priors.Norm_Prec*gmm.priors.Norm_Mu);
      % covariance posterior
      gmm.post(k).Wish_alpha=0.5*(gmm.priors.Wish_alpha+Nm(k));
      gmm.post(k).Wish_B=gmm.priors.Wish_B+Nm(k)*Cov;
      gmm.post(k).Wish_iB=inv(gmm.post(k).Wish_B);
      % weights posterior
      gmm.post(k).Dir_alpha=gmm.priors.Dir_alpha(k)+Nm(k);
    end;
    
    % assing using k-means
   case 'kmeans'
    [centre,Sigma,weight,Nm]=kminit(data,K,options);
    % Moments of parameters
    Ecentre=mean(centre');		% Expectation of means
    Covcentre=cov(centre');		% Covariance of means
    ESigma=squeeze(sum(Sigma,3))./K;	
    % Expectation of Covariances Covariance of Covariances
    Sigmatmp=reshape(Sigma,ndim*ndim,K)';
    CovSigma=reshape(cov(Sigmatmp),ndim*ndim,ndim*ndim); 
    % Initialising the posteriors
    for k=1:K,
      gmm.post(k).Norm_Mu=centre(:,k);
      gmm.post(k).Norm_Cov=Covcentre;
      Sigmatmp=squeeze(Sigma(:,:,k));
      alpha=Nm(k);
      gmm.post(k).Wish_alpha=alpha;
      gmm.post(k).Wish_B=Sigmatmp*alpha;
      gmm.post(k).Wish_iB=inv(gmm.post(k).Wish_B);
      gmm.post(k).Dir_alpha=Nm(k);
    end;
   otherwise
    error('Unknown intialisation option');
  end;
 
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [centre,Sigma,weight,Nm]=kminit(data,K,options)
% k-means for initialisation
% [centre,Sigma,weight,Nm]=kminit(data,K,options)
%
% centre      cluster means
% Sigma       cluster covariances
% weight      cluster weight
% Nm          cluster cardinality


[N,ndim]=size(data);

[mue,clv,HV] = k_means(data,K);
% get the means
centre=mue;
% compute the variances
for k=1:K,
  ndx=find(clv==k);
  Nm(k)=length(ndx);
  if Nm(k)<ndim,
    error(['Initialisation: too few samples for Covariance' ...
	   ' computation. Please Restart']);
  else
    Sigmatmp=cov(data(find(clv==k),:));
    err=testcovmatrix(Sigmatmp);
    if err.number & options.testcovmat
      %error(['Initialisation: ' err.message]);
      Sigmatmp=Sigmatmp+eye(size(Sigmatmp));
    end;
    Sigma(:,:,k)=Sigmatmp;
  end;
  % computing the weights
  weight(k)=Nm(k)/N;	% ratio of classes in total sample
end;
