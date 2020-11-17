function [pValM, pValS] = computeSignificance(dataSet1, dataSet2, noiseMean, noiseStd, initMu1, initSig1, initMu2, initSig2)

    pValM = 0.5;
    pValS = 0.5;
    return;

    dataSet1 = dataSet1(:);
    dataSet2 = dataSet2(:);
    useRejection = true;
    
    if useRejection
        N = 4000000;
        
        minMu = -2;
        maxMu = max([initMu1+5,initMu2+5]);
        
        minSig = min([initSig1*0.5,initSig2*0.5]);
        maxSig = max([initSig1*5.0,initSig2*5.0]);
        minSig = max(minSig, 0.1);
        
        mu1_samples = unifrnd(minMu, maxMu,[1,N]);
        mu2_samples = unifrnd(minMu, maxMu,[1,N]);
        sig1_samples = unifrnd(minSig, maxSig,[1,N]);
        sig2_samples = unifrnd(minSig, maxSig,[1,N]);
        
        resConv = 2e3;

        
        %% First data set
        ppm = ParforProgressbar(N, 'parpool', {'local', 8});
        done = zeros(1,N);
        logprobs1 = zeros(1,N);
        parfor sampIdx = 1:N
            %for d1Idx = 1:size(dataSet1,1)
                [logprobs_temp, norm_const] = sumLognormNormpdf(dataSet1, mu1_samples(sampIdx), sig1_samples(sampIdx), noiseMean, noiseStd,resConv);
                logprobs1(sampIdx) = sum(logprobs_temp) - log(norm_const)*size(logprobs_temp,1);
            %end
            %disp(['Done ', num2str(sampIdx)]);
            ppm.increment();

        end

        delete(ppm);

        probs1 = exp(logprobs1);
        
        remNum = 0;
        for k=1:remNum
            probs1 = probs1/max(probs1);
            valid = probs1 < 1;
            probs1 = probs1(valid);
            mu1_samples = mu1_samples(valid);
            sig1_samples = sig1_samples(valid);
        end

        probs1 = probs1/max(probs1);
        accept1 = rand(size(probs1)) < probs1;
        mu1_samples_rej = mu1_samples(accept1);
        sig1_samples_rej = sig1_samples(accept1);
        
        %% second data set
        ppm = ParforProgressbar(N, 'parpool', {'local', 8});
        done = zeros(1,N);
        logprobs2 = zeros(1,N);
        parfor sampIdx = 1:N
            %for d1Idx = 1:size(dataSet1,1)
                [logprobs_temp, norm_const] = sumLognormNormpdf(dataSet2, mu2_samples(sampIdx), sig2_samples(sampIdx), noiseMean, noiseStd,resConv);
                logprobs2(sampIdx) = sum(logprobs_temp) - log(norm_const)*size(logprobs_temp,1);
            %end
            %disp(['Done ', num2str(sampIdx)]);
            ppm.increment();

        end

        delete(ppm);

        probs2 = exp(logprobs2);
        
        remNum = 0;
        for k=1:remNum
            probs2 = probs2/max(probs2);
            valid = probs2 < 1;
            probs2 = probs2(valid);
            mu2_samples = mu2_samples(valid);
            sig2_samples = sig2_samples(valid);
        end

        probs2 = probs2/max(probs2);
        accept2 = rand(size(probs2)) < probs2;
        mu2_samples_rej = mu2_samples(accept2);
        sig2_samples_rej = sig2_samples(accept2);
        
        %%
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
            subplot(1,2,1);
            histogram(mu1_samples_rej,muEdges);
            title('mean dataset 1');

            subplot(1,2,2);
            histogram(mu2_samples_rej,muEdges);
            title('mean dataset 2');

            figure;
            subplot(1,2,1);
            histogram(sig1_samples_rej,sigEdges);
            title('std dataset 1');

            subplot(1,2,2);
            histogram(sig2_samples_rej,sigEdges);
            title('std dataset 2');
        end
        
        
    end

    useStan = false;
    if useStan
        nSamples = 100;
        
        initSig1 = initSig1;
        stan_data = struct('Nsamples1',size(dataSet1,1),'data1',dataSet1, 'noiseMu', noiseMean, 'noiseSigma', noiseStd,...
                           'mu1_lb', initMu1-5, 'mu1_ub', initMu1+5, 'sig1_lb', initSig1*0.5, 'sig1_ub', initSig1*1.5);

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
        %control.adapt_delta = 0.999;
        %control.adapt_gamma = 0.03;

        sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'init', initVals, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);%, 'control', control);
        %sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);
        %sm.is_compiled = true;
        %sm.compile();

        % subsequent calls will skip recompilation
        fit = sm.sampling('data',stan_data);

        %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
        %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'chains', 4, 'iter', 1000);
        fit.block()

        samples1 = fit.extract('permuted',true);
        
        
        subplot(1,2,1);
        histogram(samples1.mu1,200);
        
        %%
        initSig2 = initSig2;
        stan_data = struct('Nsamples1',size(dataSet2,1),'data1',dataSet2, 'noiseMu', noiseMean, 'noiseSigma', noiseStd,...
                           'mu1_lb', initMu2-5, 'mu1_ub', initMu2+5, 'sig1_lb', initSig2*0.5, 'sig1_ub', initSig2*1.5);

        
        initVals.mu1 = initMu2;
        initVals.sigma1 = initSig2;
        %initVals.mu2 = initMu2;
        %initVals.sigma2 = initSig2;


        sm2 = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'init', initVals, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);%, 'control', control);
        %sm = StanModel('model_code',model_code, 'model_name', 'significance_comp','verbose',true, 'chains', 1, 'iter', nSamples, 'file_overwrite', true);
        %sm.is_compiled = true;
        %sm.compile();

        % subsequent calls will skip recompilation
        fit2 = sm2.sampling('data',stan_data);

        %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
        %fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'chains', 4, 'iter', 1000);
        fit2.block()
        
        samples2 = fit2.extract('permuted',true);

        subplot(1,2,1);
        histogram(samples1.mu1,200);

        subplot(1,2,2);
        histogram(samples2.mu1,200);

        a = 1;
    end
end