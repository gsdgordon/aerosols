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


clipNegatives = true; %Negative counts values set to zero
%% fit FG after the event
% Integrate to get in the same window
windowSize = 20;
buffer = 20;
windowStart = -5; %Should really be 0 for all cases

sampleValid = time >= windowStart & time < (windowStart + windowSize + buffer);
sample_preInt = fg(:,sampleValid,:);
sample_int_fg_after_all = zeros(size(fg,1),size(fg,3));

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
    sample_int_fg_after_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_fg_after_all(:,1)));

sample_int_fg_after = zeros(nValid,size(sample_int_fg_after_all,2));
diameters_2 = zeros(size(sample_int_fg_after,1)+1, size(sample_int_fg_after,2));

for k = 1:size(sample_int_fg_after_all,2)
    temp = sample_int_fg_after_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_fg_after(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

if clipNegatives
    sample_int_fg_after(sample_int_fg_after<0) = 0;
end

figure;
for k=1:size(sample_int_fg_after,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_fg_after(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.01), 'mu_UB', log(50), 'sig_LB', 0.2, 'sig_UB', 50);
    A_fg_after(k) = A_t;
    mu_fg_after(k) = mu_t;
    sigma_fg_after(k) = sigma_t;
    
end

%% fit FG before the event
windowSize = 20;
buffer = 20;
windowEnd = windowStart; %Should really be 0 for all cases

sampleValid = time >= (windowStart- windowSize - buffer) & time < (windowEnd);
sample_preInt = fg(:,sampleValid,:);
sample_int_fg_before_all = zeros(size(fg,1),size(fg,3));

for k = 1:size(sample_preInt,3)
    cumdt = 0;
    cumdd = zeros(size(sample_preInt,1),1);
    for kk=1:size(sample_preInt,2)
        cumdt = cumdt + 1;
        
        if ~all(isnan(sample_preInt(:,kk,k)))
            cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);
            
            cumdt = 0;
        end
       
    end
    sample_int_fg_before_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_fg_before_all(:,1)));

sample_int_fg_before = zeros(nValid,size(sample_int_fg_before_all,2));
diameters_2 = zeros(size(sample_int_fg_before,1)+1, size(sample_int_fg_before,2));

for k = 1:size(sample_int_fg_before_all,2)
    temp = sample_int_fg_before_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_fg_before(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

if clipNegatives
    sample_int_fg_before(sample_int_fg_before<0) = 0;
end

figure;
for k=1:size(sample_int_fg_before,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_fg_before(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.01), 'mu_UB', log(50), 'sig_LB', 0.2, 'sig_UB', 50);
    A_fg_before(k) = A_t;
    mu_fg_before(k) = mu_t;
    sigma_fg_before(k) = sigma_t;
    
end

%% fit BG after the event
% Integrate to get in the same window
windowSize = 80;
buffer = 20;
windowStart = -5; %Should really be 0 for all cases

sampleValid = time >= windowStart & time < (windowStart + windowSize + buffer);
sample_preInt = bg(:,sampleValid,:);
sample_int_bg_after_all = zeros(size(bg,1),size(bg,3));

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
    sample_int_bg_after_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_bg_after_all(:,1)));

sample_int_bg_after = zeros(nValid,size(sample_int_bg_after_all,2));
diameters_2 = zeros(size(sample_int_bg_after,1)+1, size(sample_int_bg_after,2));

for k = 1:size(sample_int_bg_after_all,2)
    temp = sample_int_bg_after_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_bg_after(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

figure;
for k=1:size(sample_int_bg_after,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_bg_after(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.01), 'mu_UB', log(0.2), 'sig_LB', 0.1, 'sig_UB', 50);
    A_bg_after(k) = A_t;
    mu_bg_after(k) = mu_t;
    sigma_bg_after(k) = sigma_t;
    
end

%% fit BG before the event
windowSize = 80;
buffer = 20;
windowEnd = windowStart; %Should really be 0 for all cases

sampleValid = time >= (windowStart- windowSize - buffer) & time < (windowEnd);
sample_preInt = bg(:,sampleValid,:);
sample_int_bg_before_all = zeros(size(bg,1),size(bg,3));

for k = 1:size(sample_preInt,3)
    cumdt = 0;
    cumdd = zeros(size(sample_preInt,1),1);
    for kk=1:size(sample_preInt,2)
        cumdt = cumdt + 1;
        
        if ~all(isnan(sample_preInt(:,kk,k)))
            cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);
            
            cumdt = 0;
        end
       
    end
    sample_int_bg_before_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_bg_before_all(:,1)));

sample_int_bg_before = zeros(nValid,size(sample_int_bg_before_all,2));
diameters_2 = zeros(size(sample_int_bg_before,1)+1, size(sample_int_bg_before,2));

for k = 1:size(sample_int_bg_before_all,2)
    temp = sample_int_bg_before_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_bg_before(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

figure;
for k=1:size(sample_int_bg_before,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_bg_before(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.01), 'mu_UB', log(0.2), 'sig_LB', 0.1, 'sig_UB', 50);
    A_bg_before(k) = A_t;
    mu_bg_before(k) = mu_t;
    sigma_bg_before(k) = sigma_t;
    
end
%% Plot
subplot(3,2,1)
plot(A_fg_before,A_fg_after,'x');
hold on;
plot(A_fg_before,A_fg_before);
title('Amount of aerosol');

subplot(3,2,3)
plot(mu_fg_before,mu_fg_after,'x');
hold on;
plot(mu_fg_before,mu_fg_before);
title('mu');

subplot(3,2,5)
plot(sigma_fg_before,sigma_fg_after,'x');
hold on;
plot(sigma_fg_before,sigma_fg_before);
title('sigma');

subplot(3,2,2)
plot(A_bg_before,A_bg_after,'x');
hold on;
plot(A_bg_before,A_bg_before);
title('Amount of aerosol');

subplot(3,2,4)
plot(mu_bg_before,mu_bg_after,'x');
hold on;
plot(mu_bg_before,mu_bg_before);
title('mu');

subplot(3,2,6)
plot(sigma_bg_before,sigma_bg_after,'x');
hold on;
plot(sigma_bg_before,sigma_bg_before);
title('sigma');

%% Fits
% [A_fg_before_mu, A_fg_before_sig] = normfit(A_fg_before);
% [A_fg_after_mu, A_fg_after_sig] = normfit(A_fg_after);
% [A_bg_before_mu, A_bg_before_sig] = normfit(A_bg_before);
% [A_bg_after_mu, A_bg_after_sig] = normfit(A_bg_after);
% A_fg_before_phat = gamfit(A_fg_before);
% A_fg_after_phat = gamfit(A_fg_after);
% A_bg_before_phat = gamfit(A_bg_before);
% A_bg_after_phat = gamfit(A_bg_after);
A_fg_before_phat = gamfit(A_fg_before);
A_fg_after_phat = gamfit(A_fg_after);
A_bg_before_phat = lognfit(A_bg_before);
A_bg_after_phat = lognfit(A_bg_after);


[mu_fg_before_mu, mu_fg_before_sig] = normfit(mu_fg_before);
[mu_fg_after_mu, mu_fg_after_sig] = normfit(mu_fg_after);
[mu_bg_before_mu, mu_bg_before_sig] = normfit(mu_bg_before);
[mu_bg_after_mu, mu_bg_after_sig] = normfit(mu_bg_after);

sigma_fg_before_phat = gamfit(sigma_fg_before);
sigma_fg_after_phat = gamfit(sigma_fg_after);
sigma_bg_before_phat = gamfit(sigma_bg_before);
sigma_bg_after_phat = gamfit(sigma_bg_after);

plot_A = linspace(0,1e7,500);
plot_mu = linspace(-5,1,500);
plot_sig = linspace(0,2,500);

figure;
subplot(6,2,1);
plot(exp(plot_mu), normpdf(plot_mu,mu_fg_before_mu, mu_fg_before_sig));
title('mu fg before dist');

subplot(6,2,2);
plot(exp(plot_mu), normpdf(plot_mu,mu_fg_after_mu, mu_fg_after_sig));
title('mu fg after dist');

subplot(6,2,3);
plot(exp(plot_sig), gampdf(plot_sig,sigma_fg_before_phat(1), sigma_fg_before_phat(2)));
title('sigma fg before dist');

subplot(6,2,4);
plot(exp(plot_sig), gampdf(plot_sig,sigma_fg_after_phat(1), sigma_fg_after_phat(2)));
title('sigma fg after dist');

subplot(6,2,5);
plot((plot_A), gampdf(plot_A,A_fg_before_phat(1), A_fg_before_phat(2)));
title('A fg before dist');

subplot(6,2,6);
plot((plot_A), gampdf(plot_A,A_fg_after_phat(1), A_fg_after_phat(2)));
title('A fg after dist');

subplot(6,2,7);
plot(exp(plot_mu), normpdf(plot_mu,mu_bg_before_mu, mu_bg_before_sig));
title('mu bg before dist');

subplot(6,2,8);
plot(exp(plot_mu), normpdf(plot_mu,mu_bg_after_mu, mu_bg_after_sig));
title('mu bg after dist');

subplot(6,2,9);
plot(exp(plot_sig), gampdf(plot_sig,sigma_bg_before_phat(1), sigma_bg_before_phat(2)));
title('sigma bg before dist');

subplot(6,2,10);
plot(exp(plot_sig), gampdf(plot_sig,sigma_bg_after_phat(1), sigma_bg_after_phat(2)));
title('sigma bg after dist');

subplot(6,2,11);
plot((plot_A), lognpdf(plot_A,A_bg_before_phat(1), A_bg_before_phat(2)));
title('A bg before dist');

subplot(6,2,12);
plot((plot_A), lognpdf(plot_A,A_bg_after_phat(1), A_bg_after_phat(2)));
title('A bg after dist');




%% Now run STAN code
stan_data = struct('Nbins',size(diameters_2,1),'Nsamples',size(sample_int_fg_after,2),'data_counts',sample_int_fg_after','diameters',diameters_2');

N = size(sample_int_fg_after,2);
alpha_est = N*sum(sigma_fg_after)/(N*sum(sigma_fg_after .* log(sigma_fg_after)) - sum(log(sigma_fg_after)*sum(sigma_fg_after)));
theta_est = 1/N^2 * (N*sum(sigma_fg_after .* log(sigma_fg_after)) - sum(log(sigma_fg_after)*sum(sigma_fg_after)));
beta_est = 1/theta_est;


initVals.mu_mu = mean(mu_fg_after) - log(2);
initVals.mu_sigma = std(mu_fg_after);

phat = gamfit(sigma_fg_after);
initVals.sigma_alpha = phat(1);
initVals.sigma_beta = phat(2);

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

fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
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
