function [pValM, pValS, samples1, samples2] = computeSignificance(dataSet1, dataSet2, noiseMean, noiseStd, initMu1, initSig1, initMu2, initSig2, varargin)

    %pValM = 0.04;
    %pValS = 0.04;
    %return;
    
    p = inputParser;
    addOptional(p,'muDeltaUB',5);
    addOptional(p,'muMinIn',-2);
    addOptional(p,'sigFacU',5.0);
    addOptional(p,'sigFacL',0.5);
    addOptional(p,'inputSamples1',[]);
    addOptional(p,'inputSamples2',[]);
    
    parse(p,varargin{:});

    muDeltaUB = p.Results.muDeltaUB;
    muMinIn = p.Results.muMinIn;
    sigFacU = p.Results.sigFacU;
    sigFacL = p.Results.sigFacL;
    inputSamples1 = p.Results.inputSamples1;
    inputSamples2 = p.Results.inputSamples2;

    dataSet1 = dataSet1(:);
    dataSet2 = dataSet2(:);
    
    useRejection = true;
    
    if useRejection
        N_start = 2000000;
        
        %minMu = log(noiseMean + 1*noiseStd);
        minMu = muMinIn;
        maxMu = max([initMu1+muDeltaUB,initMu2+muDeltaUB]);
        
        minSig = min([initSig1*sigFacL,initSig2*sigFacL]);
        maxSig = max([initSig1*sigFacU,initSig2*sigFacU]);
        minSig = max(minSig, 0.1);
        maxSig = max(maxSig,minSig+0.1);
        
        resConv = 2e3;
        
        minSamples = 500;
        
        %% First data set
        if isempty(inputSamples1)
            nSamples = 0;
            mu1_samples_rej = [];
            sig1_samples_rej = [];
            N_curr = N_start;
            while (nSamples < minSamples)

                mu1_samples = unifrnd(minMu, maxMu,[1,N_curr]);
                mu1_samples(end) = initMu1;
                dmu1 = exp(mu1_samples);
                %mu1_samples = log(unifrnd(exp(minMu), exp(maxMu),[1,N_curr]));
                sig1_samples = unifrnd(minSig, maxSig,[1,N_curr]);
                sig1_samples(end) = initSig1;

                ppm = ParforProgressbar(N_curr, 'parpool', {'local', 8});
                logprobs1 = zeros(1,N_curr);
                parfor sampIdx = 1:N_curr
                    %for d1Idx = 1:size(dataSet1,1)
                        [logprobs_temp, norm_const] = sumLognormNormpdf(dataSet1, mu1_samples(sampIdx), sig1_samples(sampIdx), noiseMean, noiseStd,resConv, exp(log(max(dataSet1))+3*initSig1));
                        logprobs1(sampIdx) = sum(logprobs_temp) - log(norm_const)*size(logprobs_temp,1);
                    %end
                    %disp(['Done ', num2str(sampIdx)]);
                    ppm.increment();

                end

                delete(ppm);

                logprobs1 = logprobs1 + log(dmu1);
                validProbs = ~isinf(logprobs1) & ~isnan(logprobs1);% & (logprobs1 < 0);
                logprobs1 = logprobs1(validProbs);

                probNorm_log = max(logprobs1);

                logprobs1 = logprobs1 - probNorm_log;
                probs1 = exp(logprobs1);
                accept1 = rand(size(probs1)) < probs1;
                mu1_samples_rej = mu1_samples(accept1);
                sig1_samples_rej = sig1_samples(accept1);

                nSamples = size(sig1_samples_rej,2);

                incFac = minSamples/nSamples;
                N_curr = ceil(1.05*(incFac) * N_curr);
                
                %N_curr = min(N_curr, N_start*10);
                
                %minMu_curr = min(mu1_samples_rej);
                %maxMu_curr = max(mu1_samples_rej);
                
                
            end
        else
            mu1_samples_rej = inputSamples1(:,1);
            sig1_samples_rej = inputSamples1(:,2);
        end
        
        %samples1 = [mu1_samples_rej(:), sig1_samples_rej(:)];
        
        %% second data set
        if isempty(inputSamples2)
            nSamples = 0;
            mu2_samples_rej = [];
            sig2_samples_rej = [];
            N_curr = N_start;
            while (nSamples < minSamples)
                mu2_samples = unifrnd(minMu, maxMu,[1,N_curr]);
                mu2_samples(end) = initMu2;
                dmu2 = exp(mu2_samples);
                %mu2_samples = log(unifrnd(exp(minMu), exp(maxMu),[1,N_curr]));
                sig2_samples = unifrnd(minSig, maxSig,[1,N_curr]);
                sig2_samples(end) = initSig2;

                ppm = ParforProgressbar(N_curr, 'parpool', {'local', 8});
                logprobs2 = zeros(1,N_curr);
                parfor sampIdx = 1:N_curr
                    %for d1Idx = 1:size(dataSet1,1)
                        [logprobs_temp, norm_const] = sumLognormNormpdf(dataSet2, mu2_samples(sampIdx), sig2_samples(sampIdx), noiseMean, noiseStd,resConv, exp(log(max(dataSet2))+3*initSig2));
                        logprobs2(sampIdx) = sum(logprobs_temp) - log(norm_const)*size(logprobs_temp,1);
                    %end
                    %disp(['Done ', num2str(sampIdx)]);
                    ppm.increment();

                end

                delete(ppm);

                logprobs2 = logprobs2 + log(dmu2);
                validProbs = ~isinf(logprobs2) & ~isnan(logprobs2);% & (logprobs2 < 0);
                logprobs2 = logprobs2(validProbs);
                probNorm_log = max(logprobs2);

                logprobs2 = logprobs2 - probNorm_log;
                probs2 = exp(logprobs2);
                accept2 = rand(size(probs2)) < probs2;
                mu2_samples_rej = mu2_samples(accept2);
                sig2_samples_rej = sig2_samples(accept2);

                nSamples = size(sig2_samples_rej,2);
                incFac = minSamples/nSamples;
                N_curr = ceil(1.05*(incFac) * N_curr);
                
                %N_curr = min(N_curr, N_start*10);
            end
        else
            mu2_samples_rej = inputSamples2(:,1);
            sig2_samples_rej = inputSamples2(:,2);
        end
        
        %samples2 = [mu2_samples_rej(:), sig2_samples_rej(:)];
        
        %%
        
        
    end

    useStan = false;
    if useStan
        nSamples = 1000;
        
        minMu = muMinIn;
        maxMu = max([initMu1+muDeltaUB,initMu2+muDeltaUB]);
        
        minSig = min([initSig1*sigFacL,initSig2*sigFacL]);
        maxSig = max([initSig1*sigFacU,initSig2*sigFacU]);
        minSig = max(minSig, 0.1);
        maxSig = max(maxSig,minSig+0.1);
        
        initSig1 = initSig1;
        stan_data = struct('Nsamples1',size(dataSet1,1),'data1',dataSet1, 'noiseMu', noiseMean, 'noiseSigma', noiseStd,...
                           'mu1_lb', minMu, 'mu1_ub', maxMu, 'sig1_lb', minSig, 'sig1_ub', maxSig);

        initVals.mu1 = initMu1;
        initVals.sigma1 = initSig1;
        %initVals.mu2 = initMu2;
        %initVals.sigma2 = initSig2;

        fileID = fopen('significance_comp.stan');
        model_code = textscan(fileID,'%s');
        model_code = model_code{:};
        fclose(fileID);

        %control.stepsize = 3;
        %control.stepsize_jitter = 1;
        control.adapt_delta = 0.999;
        %control.adapt_gamma = 0.03;

        nChains = 8;
        sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'init', initVals, 'chains', nChains, 'iter', nSamples/nChains, 'file_overwrite', true, 'warmup', 5*nSamples/nChains, 'control', control);
        %sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);
        %sm.is_compiled = true;
        %sm.compile();

        if isempty(inputSamples1)
            % subsequent calls will skip recompilation
            fit = sm.sampling('data',stan_data);

            %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
            %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'chains', 4, 'iter', 1000);
            fit.block()

            samples1 = fit.extract('permuted',true);
            
            %subplot(1,2,1);
            %histogram(samples1.mu1,200);
            
            
            mu1_samples_rej = samples1.mu1;
            sig1_samples_rej = samples1.sigma1;
        else
            mu1_samples_rej = inputSamples1(:,1);
            sig1_samples_rej = inputSamples1(:,2);
        end
        

        
        %%
        initSig2 = initSig2;
        stan_data = struct('Nsamples1',size(dataSet2,1),'data1',dataSet2, 'noiseMu', noiseMean, 'noiseSigma', noiseStd,...
                           'mu1_lb', initMu2-5, 'mu1_ub', initMu2+5, 'sig1_lb', initSig2*0.5, 'sig1_ub', initSig2*1.5);

        
        initVals.mu1 = initMu2;
        initVals.sigma1 = initSig2;
        %initVals.mu2 = initMu2;
        %initVals.sigma2 = initSig2;

        sm2 = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'init', initVals, 'chains', nChains, 'iter', nSamples/nChains, 'file_overwrite', true, 'warmup', 5*nSamples/nChains, 'control', control);
        %sm2 = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'init', initVals, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);%, 'control', control);
        %sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);
        %sm.is_compiled = true;
        %sm.compile();

        

        if isempty(inputSamples2)
            % subsequent calls will skip recompilation
            fit2 = sm2.sampling('data',stan_data);

            %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
            %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'chains', 4, 'iter', 1000);
            fit2.block()

            samples2 = fit2.extract('permuted',true);

            mu2_samples_rej = samples2.mu1;
            sig2_samples_rej = samples2.sigma1;
        else
            mu2_samples_rej = inputSamples2(:,1);
            sig2_samples_rej = inputSamples2(:,2);
        end

        %subplot(1,2,1);
        %histogram(samples1.mu1,200);

        %subplot(1,2,2);
        %histogram(samples2.mu1,200);

        a = 1;
    end
    
    samples1 = [mu1_samples_rej(:), sig1_samples_rej(:)];
    samples2 = [mu2_samples_rej(:), sig2_samples_rej(:)];
    
    minVal = min([min(mu1_samples_rej), min(mu2_samples_rej)]);
    maxVal = max([max(mu1_samples_rej), max(mu2_samples_rej)]);
    muEdges = linspace(minVal-1e-20, maxVal+1e-20,50+1);

    if size(mu1_samples_rej,2) > 10000
        mu1_samples_rej = mu1_samples_rej(1, randperm(size(mu1_samples_rej,2),10000));
        sig1_samples_rej = sig1_samples_rej(1, randperm(size(mu1_samples_rej,2),10000));
    end

    if size(mu2_samples_rej,2) > 10000
        mu2_samples_rej = mu2_samples_rej(1, randperm(size(mu2_samples_rej,2),10000));
        sig2_samples_rej = sig2_samples_rej(1, randperm(size(mu2_samples_rej,2),10000));
    end

    [M1, M2] = meshgrid(mu1_samples_rej, mu2_samples_rej);
    countMatM = M2 > M1;
    pValM = nnz(countMatM)/numel(countMatM);

    minVal = min([min(sig1_samples_rej), min(sig2_samples_rej)]);
    maxVal = max([max(sig1_samples_rej), max(sig2_samples_rej)]);
    sigEdges = linspace(minVal-1e-20, maxVal+1e-20,50+1);

    [S1, S2] = meshgrid(sig1_samples_rej, sig2_samples_rej);
    countMatS = S2 > S1;
    pValS = nnz(countMatS)/numel(countMatS);


    plotHists = false;

    if plotHists
        figure;
        subplot(2,2,1);
        histogram(mu1_samples_rej,muEdges);
        title('mean dataset 1');

        subplot(2,2,2);
        histogram(mu2_samples_rej,muEdges);
        title('mean dataset 2');

        subplot(2,2,3);
        histogram(sig1_samples_rej,sigEdges);
        title('std dataset 1');

        subplot(2,2,4);
        histogram(sig2_samples_rej,sigEdges);
        title('std dataset 2');

        a = 1;
    end
end