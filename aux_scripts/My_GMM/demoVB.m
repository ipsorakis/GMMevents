clf; clc;
fprintf('Loading the data file\n');
load gmmEg;

% fprintf('\nPlot the data'); 
plotset(x, trueLabels,4);
fprintf('\npress any key to continue'); pause;clc;

K=100;
fprintf('Now run the variational Bayes algorithm with K=%d components - so very over-complex model\n\n', K);
[gmm,fe] = gmmvar(x,K);
fprintf('\nShow the posterior expectations of mixings, p(k)\n');
clf;
bar(mean(gmm.pjgx));
title('Only 3 components?','fontsize',15);
fprintf('\npress any key to continue'); pause; clc;

fprintf('Now draw samples from the distributions over the means and covariances \n');
fprintf('The means have m ~ N(m_mu,m_cov)  and the covariances C ~ Wi(cov_B/cov_alpha , N-1)\n')
fprintf('Where N is the count hyper-parameter from the Dirichlet pdf\n');
fprintf('\npress any key to continue'); pause; clc;

clf;
plotset(x,trueLabels,4); hold on;
for k=1:gmm.K,
    Nk = floor(gmm.post(k).Dir_alpha - 1);
    if (Nk > 2), % the wishart is singular otherwise
        for n=1:10,
            mk = mvnrnd(gmm.post(k).Norm_Mu,gmm.post(k).Norm_Cov)'; % draw a mean sample
            Ck = wishrnd(gmm.post(k).Wish_B./(gmm.post(k).Wish_alpha*Nk),Nk); % draw a covariance sample
            plotellipse(mk,Ck,'c');
            plot(mk(1),mk(2),'c.','markersize',20);
        end;
        mk = gmm.post(k).Norm_Mu; % E{p(mu)}
        Ck = gmm.post(k).Wish_B./gmm.post(k).Wish_alpha; % E{p(C)}
        plotellipse(mk,Ck,'m');
        plot(mk(1),mk(2),'m.','markersize',20);
    end;
end;
title({'Cyan - means and covariance posterior samples.', 'Magenta - expectation over posteriors'},'fontsize',15);
hold off;