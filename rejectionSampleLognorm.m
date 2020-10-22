function samples = rejectionSampleLognorm(counts, diameters, mu_mu_est, mu_sigma_est, sigma_alpha_est, sigma_beta_est, mu_ml, sigma_ml)


    radii = diameters/2;
    log_radii = log(radii);

    log_radii_av = (log_radii(1:end-1)+log_radii(2:end))/2;
    log_bin_sizes = log(log_radii(2:end)) - log(log_radii(1:end-1));

    nStartSamples = 10000;
    
    %mu_samples = normrnd(mu_mu_est, mu_sigma_est, nStartSamples,1);
    mu_samples = unifrnd(mu_mu_est-10*mu_sigma_est,mu_mu_est+10*mu_sigma_est,nStartSamples,1);
    %sig_samples = gamrnd(sigma_alpha_est, sigma_beta_est, nStartSamples,1);
    sig_samples = unifrnd(0,sigma_alpha_est/sigma_beta_est + 6*sqrt(sigma_alpha_est)/sigma_beta_est,nStartSamples,1);
    
    nIntSteps = 4;
    intFun = @(y) trapz(normpdf(linspace(y(1),y(2),nIntSteps),y(3),y(4)))*(y(2)-y(1))/(nIntSteps-1); % For particle count
    
    nData = size(counts,2);
    overallProb = zeros(nStartSamples,1);
    
    for k=1:nData
    	currentData = counts(:,k);
        prob_ml(k) = normApproxMNpdf([mu_ml(k) - log(2),sigma_ml(k)],log_radii(:,k)',currentData',intFun);
    end
    
    w = warning;
    warning('off','all');
    for m=1:nStartSamples
        probs = zeros(nData,1);
        for k=1:nData
            currentData = counts(:,k);
            probs(k,1) = normApproxMNpdf([mu_samples(m,1),sig_samples(m)],log_radii(:,k)',currentData',intFun);
            
            probs(k,1) = probs(k,1) - prob_ml(k);
            %probs(k,1) = probs(k,1) - max(prob_ml);
        end
        
        probs(probs > 0) = -inf;
        %probs(probs < -1e20) = -inf;
        overallProb(m) = max(probs);

        
    end
    warning(w);
    

    linProb = exp(overallProb);
    accept = logical(binornd(1,linProb));
    
    mu_samples_sampled = mu_samples(accept);

    a = 1;
    
end