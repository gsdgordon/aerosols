% File to process multipe aerotrak readings for some particular event
% lognormal distribution and creat a video
% Author: George Gordon
% Date: 30/09/2020


% Clear any old data
clc;
clear variables;
close all;

%addpath('./MatlabProcessManager-master');
%addpath('./MatlabStan-2.15.1.0');
%addpath('/home/george/Apps/cmdstan');

% Load data
opts = detectImportOptions(['extubation.csv']);
opts.DataLines = [1, 5];
opts.RowNamesColumn = 1;
opts.VariableNamesLine = 0;
opts.PreserveVariableNames = true;
opts.VariableTypes{1} = 'char';
opts.VariableTypes{2} = 'char';
headers = readtable(['extubation.csv'],opts, 'ReadRowNames', true, 'ReadVariableNames', false);
headers = headers(1:5,2);

startTime = str2double(headers(3,1).Var2{1});
endTime = str2double(headers(4,1).Var2{1});
timeStep = str2double(headers(5,1).Var2{1});

T = readtable(['extubation.csv'],'ReadVariableNames', true, 'HeaderLines',6);

indices = unique(T.Index);
nItems = size(indices,1);

maxDiameter = 25; % Should load this from datasheet?
diameters = T(T.Index == indices(1),26);
diameters = table2array(diameters);
diameters = [diameters; maxDiameter];

time = startTime:timeStep:endTime;
nTimes = (endTime - startTime)/timeStep + 1;

validRows = T.Index == 1; % FIX should loop for other indices and check
tempData = table2array(T(validRows,27:27+nTimes-1));
bg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
fg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
avSampleTimes = zeros(size(indices,1),1);

for currentIdx = indices'
    % Step to check validity of data? Make sure volume/m^3 is correct
    
    validRows = T.Index == currentIdx;

    data = table2array(T(validRows,27:27+nTimes-1));
    
    % Find valid diameters
    validDiams = ~all(isnan(data),2);
    
    currentDiams_temp = diameters([validDiams; true]);
    currentDiams_av_temp = (currentDiams_temp(1:end-1)+currentDiams_temp(2:end))/2; % Mean diameter in each bin

    % Convert counts to volumes
    vols_temp = 4/3*pi*(currentDiams_av_temp/2).^3 * (1e-6)^3;
    bin_sizes_temp = currentDiams_temp(2:end) - currentDiams_temp(1:end-1);
    log_bin_sizes_temp = log(currentDiams_temp(2:end)) - log(currentDiams_temp(1:end-1));
    
    
    currentDiams = nan(size(data,1),1);
    currentDiams(validDiams) = currentDiams_temp(1:end-1);
    currentDiams_av = nan(size(data,1),1);
    currentDiams_av(validDiams) = currentDiams_av_temp;
    currentVols = nan(size(data,1),1);
    currentVols(validDiams) = vols_temp;
    currentBinSizes = nan(size(data,1),1);
    currentBinSizes(validDiams) = bin_sizes_temp;
    currentLogBinSizes = nan(size(data,1),1);
    currentLogBinSizes(validDiams) = log_bin_sizes_temp;
    
    
    data_v = data .* repmat(currentVols,1, size(data,2));

    % Densities so that a probability density approach can be used
    data_v_density = data_v ./ repmat(currentLogBinSizes,1, size(data,2)); %Try using log binsizes
    data_density = data ./ repmat(currentLogBinSizes,1, size(data,2));
    
    tempValid = ~isnan(data);
    tempValid = nansum(tempValid,1) > 0;
    avSampleTime = median(diff(time(tempValid)));
    avSampleTimes(currentIdx) = avSampleTime;
    [bg_current, fg_current] = splitBGFG(data, avSampleTime, tempValid);
    
    nSizes = size(data,1);
    tColor = lines(nSizes);

    for k=1:nSizes

        subplot(nSizes,1,k);
        currentValid = ~isnan(data(k,:));
        plot(time(currentValid),data(k,currentValid),'Color',tColor(k,:));
        hold on;
        plot(time(currentValid),bg_current(k,currentValid),'Color','black','LineStyle',':', 'LineWidth',1);

        title(['Diameter: ', num2str(diameters(k)), '\mum']);
        ylabel('#/m^3');
        xlabel('time')

        hold on;
    end

    bg(:,:,currentIdx) = bg_current;
    fg(:,:,currentIdx) = fg_current;
end

% Integrate to get in the same window
windowStart = -5; %Should really be 0 for all cases
windowSize = 20;
buffer = 20;

% fit FG
sampleValid = time >= windowStart & time < (windowStart + windowSize + buffer);
sample_preInt = fg(:,sampleValid,:);
sample_int = zeros(size(fg,1),size(fg,3));

for k = 1:size(sample_preInt,3)
    cumdt = 0;
    cumdd = zeros(size(sample_preInt,1),1);
    for kk=2:size(sample_preInt,2) % start from 2 to exclude the first point as data is cumulative after the event
        cumdt = cumdt + 1;
        
        if ~all(isnan(sample_preInt(:,kk,k)))
            cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);
            
            cumdt = 0;
        end
       
    end
    sample_int(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int(:,1)));

sample_int_2 = zeros(nValid,size(sample_int,2));
diameters_2 = zeros(size(sample_int_2,1)+1, size(sample_int_2,2));

for k = 1:size(sample_int,2)
    temp = sample_int(:,k);
    valid = ~isnan(temp);
    
    sample_int_2(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

%% Now run STAN code
stan_data = struct('Nbins',size(diameters_2,1),'Nsamples',size(sample_int_2,2),'data_counts',sample_int_2','diameters',diameters_2');

currentLogBinSizes = log(diameters_2(2:end,1)) - log(diameters_2(1:end-1,1));
[A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,1).', (sample_int(2:end,1)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.01), 'mu_UB', log(50), 'sig_LB', 0.1, 'sig_UB', 50);
%initVals.mu = mu_t - log(2);
%initVals.sigma = sigma_t;
initVals.mu_mu = mu_t - log(2);
initVals.mu_sigma = 0.01;
scale = 1;
initVals.sigma_alpha = 2;
initVals.sigma_beta = initVals.sigma_alpha/sigma_t;

figure;
subplot(2,1,1);
plot_mu = linspace(-5,5,500);
plot_mu_mu = normpdf(plot_mu,initVals.mu_mu,initVals.mu_sigma);
plot(plot_mu, plot_mu_mu);
title('distribtion of mu');

subplot(2,1,2);
plot_sig = linspace(0,5,500);
plot_sig_sig = gampdf(plot_sig,initVals.sigma_alpha,initVals.sigma_beta);
plot(plot_sig, plot_sig_sig);
title('distribtion of sigma');

fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 8, 'iter', 125000);
%fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'chains', 4, 'iter', 1000);
fit.block()

samples = fit.extract('permuted',true);

figure;
subplot(3,2,1);
hist(samples.mu,100);
title('mu');

subplot(3,2,2);
hist(samples.sigma,100);
title('sigma');

subplot(3,2,3);
hist(samples.mu_mu,100);
title('mu mu');

subplot(3,2,4);
hist(samples.mu_sigma,100);
title('mu sigma');

subplot(3,2,5);
hist(samples.sigma_alpha,100);
title('sigma alpha');

subplot(3,2,6);
hist(samples.sigma_beta,100);
title('sigma beta');

print(fit);

a = 1;
