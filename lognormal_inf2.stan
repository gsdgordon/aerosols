functions {
    real my_mvn_lpdf(vector x, vector mu_mvn, matrix sig_mvn) {
        vector[num_elements(x)] eigs;
        real det_est;
        real prob;
        real normAttempt;
        matrix[rows(sig_mvn), cols(sig_mvn)] sig_mvn_inv;
        matrix[rows(sig_mvn), cols(sig_mvn)] sig_mvn_sym;

        sig_mvn_inv = inverse(sig_mvn * sig_mvn')*sig_mvn';
        prob = -1*0.5*(x - mu_mvn)'*sig_mvn_inv*(x-mu_mvn);

        sig_mvn_sym = 0.5*(sig_mvn + sig_mvn');
        eigs = eigenvalues_sym(sig_mvn_sym);


        det_est = 1.0;
        for (n in 1:num_elements(x)) {
            if (eigs[n]  < 10^-9 * max(eigs)) {
                eigs[n] = 1;
            }
            det_est = det_est*eigs[n];
        }

        
        normAttempt = 0.5*log(abs(det_est)*(2*pi())^num_elements(x));


        prob = prob - normAttempt;

        return prob;
    
    }
}


data {
    int Nbins;
    int Nsamples;
    
    real data_density[Nsamples,Nbins-1];
    real diameters[Nsamples,Nbins];
}

transformed data {
    real newCount = 1000;
    real totalCount[Nsamples];
    matrix[Nsamples,Nbins-1] data_counts_norm;
    real logradii[Nsamples,Nbins];

    for (m in 1:Nsamples) {
        totalCount[m] = 0;
        for (n in 1:(Nbins-1)) {
            totalCount[m] += data_counts[m,n];
        }
    }

    for (m in 1:Nsamples) {
        for (n in 1:(Nbins-1)) {
            data_counts_norm[m,n] = data_counts[m,n]/totalCount[m]*newCount;
        }
    }

    for (m in 1:Nsamples) {
        for (n in 1:Nbins) {
            logradii[m, n] = log(diameters[m,n]/2);
        }
    }
}


parameters {
    //real<lower = -2, upper = 50> mu;
    //real<lower = 0.001, upper = 50> sigma;

    
    real<lower = -2, upper = 50> mu_mu;
    real<lower = 0.001, upper = 50> mu_sigma;

    real<lower = 0> sigma_alpha;
    real<lower = 0> sigma_beta;

    real mu;
    real<lower = 0> sigma;

    real<lower = 0> errSig;
}

transformed parameters {

}

model {
    vector[Nbins-1] normprobs;
    real probNorm;
    matrix[Nbins-1,Nbins-1] P;
    matrix[Nbins-1,Nbins-1] M;
    vector[Nbins-1] mu_mvn;
    matrix[Nbins-1,Nbins-1] sig_mvn;
    matrix[Nbins-1,Nbins-1] sig_mvn_inv;

    real errmean;

    real temp;
    

    for (m in 1:Nsamples) {
        //print("DEBUG: mu_mu ", mu_mu, ", mu_sigma: ", mu_sigma);
        mu ~ normal(mu_mu, mu_sigma);

        //print("DEBUG: sigma_alpha ", sigma_alpha, ", sigma_beta: ", sigma_beta);
        sigma ~ gamma(sigma_alpha, sigma_beta);

        //print("DEBUG: mu ", mu, ", sigma: ", sigma);
        
        probNorm = 0;
        for (n in 1:(Nbins-1)) {
            errmean = normal_lcdf(logradii[m,n+1] | mu, sigma) - normal_lcdf(logradii[m, n] | mu, sigma);
            data_density[m,n] ~ normal(errmean, errsig);
        }
    }
}
