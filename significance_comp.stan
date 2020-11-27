data {
    int Nsamples1;
    //int Nsamples2;
    real sig1_lb;
    real sig1_ub;
    real mu1_lb;
    real mu1_ub;

    real data1[Nsamples1];
    //real data2[Nsamples2];

    real noiseMu;
    real noiseSigma;
}

parameters {
    real<lower = 0, upper = 1> mu1_raw;
    real<lower = 0, upper = 1> sigma1_raw;

    //real<lower = 0, upper = 100> mu2;
    //real<lower = 0, upper = 50> sigma2;

    real lognpart1_raw;
    //real lognpart1;
    real npart1_raw;
    //real lognpart2;
    //real npart2;
}

transformed parameters {
    real mu1;
    real sigma1;

    mu1 = mu1_raw * (14 - (-2)) + (-2);
    sigma1 = sigma1_raw*(10-0.1) + 0.1;
}

model {
    real lognpart1;
    real npart1;
    //sigma1 ~ uniform(sig1_lb, sig1_ub);
    //mu1 ~ uniform(mu1_lb, mu1_ub);

    for (m in 1:Nsamples1) {
        //print("DEBUG: mu1 ", mu1, ", sigma1: ", sigma1);
        //print("DEBUG: mu1 ", mu1, ", sigma1: ", sigma1);
        //lognpart1 ~ lognormal(mu1, sigma1);
        lognpart1_raw ~ normal(0,1);
        lognpart1 = exp(lognpart1_raw*sigma1 + mu1);
        
        npart1_raw ~ normal(0, 1);
        npart1 = npart1_raw*noiseSigma + noiseMu;
        data1[m] ~ normal(lognpart1 + npart1,1);
    }

    //for (m in 1:Nsamples2) {
        //print("DEBUG: mu_mu ", mu_mu, ", mu_sigma: ", mu_sigma);
    //    lognpart2 ~ lognormal(mu2, sigma2);
    //    npart2 ~ normal(noiseMu, noiseSigma);
    //    data2[m] ~ normal(lognpart2 + npart2,0.1);
    //}
}
