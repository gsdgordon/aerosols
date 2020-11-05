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

%folder = '../../StudyData';
folder = 'C:\Users\george\OneDrive - The University of Nottingham\SAVE\csv0909_1310';
file = 'UpperGI_cough.csv';
filepath = fullfile(folder,file);

% Load data
opts = detectImportOptions(filepath);

isFullData = true;
if isFullData
    opts.DataLines = [1, 5];
    opts.RowNamesColumn = 1;
    opts.VariableNamesLine = 0;
    opts.PreserveVariableNames = true;
    opts.VariableTypes{1} = 'char';
    opts.VariableTypes{2} = 'char';
    headers = readtable(filepath,opts, 'ReadRowNames', true, 'ReadVariableNames', false);
    headers = headers(1:5,2);  

    startTime = str2double(headers(3,1).Var2{1});
    endTime = str2double(headers(4,1).Var2{1});
    timeStep = str2double(headers(5,1).Var2{1});
    
    T = readtable(filepath,'ReadVariableNames', true, 'HeaderLines',5);
    
    eventTimes = table2array(unique(T(:,4)));
    indices = 1:size(eventTimes,1);
    indices = indices';
    diameters_full = T.ParticleBin_um_;
    diameters = diameters_full(table2array(T(:,4)) == eventTimes(1));
    %diameters = table2array(diameters);
    dataStartCol = 25; %% FIX should determine this from variable names
else
    T = readtable(filepath,'ReadVariableNames', true, 'HeaderLines',0);
    
    indices = 1:(size(T,1)/7);
    indices = indices';
    indices_tab = kron(indices, ones(7,1));
    T = addvars(T,indices_tab,'Before','x_180', 'NewVariableNames','Index');
    diameters = [0.3; 0.5; 0.7; 1.0; 3.0; 5.0; 10.0];
    
    startTime = -180;
    timeStep = 1;
    endTime = 180;
    
    dataStartCol = 2;
end

%% Correct for effect of tube
useTubeCorrection = true;

if useTubeCorrection
    tubeCorrection_tab = readtable('C:\Users\george\OneDrive - The University of Nottingham\SAVE\TubeCalibration\TubeBendCorrection.csv');
    %tubeCorrection_tab = readtable('/home/george/Desktop/TubeBendCorrection.csv');
    tubeCorrection = table2array(tubeCorrection_tab);

    for k=1:size(diameters,1)
        correctionIdx = find(diameters(k) == tubeCorrection(:,1));

        correctionVal(k) = tubeCorrection(correctionIdx,2);
    end
else
    correctionVal = zeros(1,nSizes);
end

%%


nItems = size(indices,1);

maxDiameter = 25; % Should load this from datasheet?
diameters = [diameters; maxDiameter];

time = startTime:timeStep:endTime;
nTimes = (endTime - startTime)/timeStep + 1;

validRows = table2array(T(:,4)) == eventTimes(1); % FIX should loop for other indices and check
tempData = table2array(T(validRows,dataStartCol:dataStartCol+nTimes-1));
bg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
fg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
avSampleTimes = zeros(size(indices,1),1);

for currentIdx = indices'
    % Step to check validity of data? Make sure volume/m^3 is correct
    
    validRows = table2array(T(:,4)) == eventTimes(currentIdx);

    data = table2array(T(validRows,dataStartCol:dataStartCol+nTimes-1));
    
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
    
    bg_current_v = bg_current .* repmat(currentVols,1, size(data,2));
    fg_current_v = fg_current .* repmat(currentVols,1, size(data,2));
    
    nSizes = size(data,1);
    tColor = lines(nSizes);

    for k=1:nSizes+2

        if (k <= nSizes)
            subplot(nSizes+2,3,3*(k-1)+1);
            currentValid = ~isnan(data(k,:));
            plot(time(currentValid),data(k,currentValid)./(1-correctionVal(k)),'Color',tColor(k,:));
            %hold on;
            %plot(time(currentValid),bg_current(k,currentValid)./(1-correctionVal(k)),'Color','black','LineStyle',':', 'LineWidth',1);

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            ylim_curr = ylim;
            if currentIdx == indices(end)
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);
            end

            if currentIdx == indices(1)
                xline(0,'k:');
            end
            hold on;

            subplot(nSizes+2,3,3*(k-1)+2);
            currentValid = ~isnan(data(k,:));
            plot(time(currentValid),bg_current(k,currentValid)./(1-correctionVal(k)),'Color',tColor(k,:),'LineStyle',':', 'LineWidth',1);

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            ylim(ylim_curr);
            if currentIdx == indices(1)
                xline(0,'k:');
            end
            hold on;

            subplot(nSizes+2,3,3*(k-1)+3);
            currentValid = ~isnan(data(k,:));
            plot(time(currentValid),fg_current(k,currentValid)./(1-correctionVal(k)),'Color',tColor(k,:));

            title(['Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time');
            ylim(ylim_curr);
            if currentIdx == indices(1)
                xline(0,'k:');
            end
            hold on;

            %pause(0.1);
        else
            if k== nSizes+1
                subplot(nSizes+2,3,3*(k-1)+1);
                currentValid = ~all(isnan(data),1);
                temp = data(:,currentValid);
                plot(time(currentValid),nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim_curr = ylim;
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;

                subplot(nSizes+2,3,3*(k-1)+2);
                temp = bg_current(:,currentValid);
                plot(time(currentValid),nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);
                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;

                subplot(nSizes+2,3,3*(k-1)+3);
                temp = fg_current(:,currentValid);
                plot(time(currentValid),nansum(fg_current(:,currentValid)./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);
                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;
            elseif k == nSizes+2
                subplot(nSizes+2,3,3*(k-1)+1);
                currentValid = ~all(isnan(data),1);
                temp = data_v(:,currentValid);
                plot(time(currentValid),nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim_curr = ylim;
                if currentIdx == indices(end)
                    ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                    ylim(ylim_curr);
                end

                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;

                subplot(nSizes+2,3,3*(k-1)+2);
                temp = bg_current_v(:,currentValid);
                plot(time(currentValid),nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim(ylim_curr);
                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;

                subplot(nSizes+2,3,3*(k-1)+3);
                temp = fg_current_v(:,currentValid);
                plot(time(currentValid),nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1),'k');

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim(ylim_curr);
                if currentIdx == indices(1)
                    xline(0,'k:');
                end
                hold on;
                
                if any(nansum(temp./repmat(1-correctionVal',1,size(temp,2)),1) > 2e-10)
                    a = 1;
                end
            end
            
        end
    end

    bg(:,:,currentIdx) = bg_current;
    fg(:,:,currentIdx) = fg_current;
    raw(:,:,currentIdx) = data;
end

disp('here');

clipNegatives = true; %Negative counts values set to zero


%% Raw difference
calcRawDiff = true;

if calcRawDiff
    rawdiffwindowSize_after = 100;
    buffer = 20;
    windowStart = 0; %Should really be 0 for all cases

    sampleValid = time > windowStart & time <= (windowStart + rawdiffwindowSize_after + buffer);
    sample_preInt = raw(:,sampleValid,:);
    sample_int_raw_after_all = zeros(size(raw,1),size(raw,3));

    for k = 1:size(sample_preInt,3)
        cumdt = 0;
        cumdd = zeros(size(sample_preInt,1),1);
        for kk=1:size(sample_preInt,2)
            
            if (kk <= rawdiffwindowSize_after)
                cumdt = cumdt + 1;
            end

            if ~all(isnan(sample_preInt(:,kk,k)))
                cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);

                cumdt = 0;
                
                if kk >= rawdiffwindowSize_after
                    break;
                end
            end

        end
        
        cumdd = cumdd/rawdiffwindowSize_after;
        sample_int_raw_after_all(:,k) = cumdd;
    end

    windowEnd = windowStart;
    rawdiffwindowSize_before = 50;
    windowStart = windowEnd - rawdiffwindowSize_before;
    
    sampleValid = time > windowStart & time <= (windowStart + rawdiffwindowSize_before + buffer);
    sample_preInt = raw(:,sampleValid,:);
    sample_int_raw_before_all = zeros(size(raw,1),size(raw,3));

    for k = 1:size(sample_preInt,3)
        cumdt = 0;
        cumdd = zeros(size(sample_preInt,1),1);
        for kk=1:size(sample_preInt,2)
            
            if (kk <= rawdiffwindowSize_before)
                cumdt = cumdt + 1;
            end

            if ~all(isnan(sample_preInt(:,kk,k)))
                cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);

                cumdt = 0;
                
                if kk >= rawdiffwindowSize_before
                    break;
                end
            end

        end
        
        cumdd = cumdd/rawdiffwindowSize_before;
        sample_int_raw_before_all(:,k) = cumdd;

    end

    sample_int_raw_diff_all = sample_int_raw_after_all - sample_int_raw_before_all;

    % Diameter array per each data set
    nValid = nnz(~isnan(sample_int_raw_diff_all(:,end))); %FIX Should be more robust than this

    sample_int_raw_diff = zeros(nValid,size(sample_int_raw_diff_all,2));
    diameters_2 = zeros(size(sample_int_raw_diff,1)+1, size(sample_int_raw_diff,2));

    for k = 1:size(sample_int_raw_diff_all,2)
        temp = sample_int_raw_diff_all(:,k);
        valid = ~isnan(temp);

        sample_int_raw_diff(:,k) = temp(valid);
        diameters_2(:,k) = diameters([valid; true]);
    end


    figure;
    for k=1:size(sample_int_raw_diff,2)
        %figure;
        currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
        currentData = sample_int_raw_diff(1:end,k);
        currentData_raw = currentData;
        currentDiameters = diameters_2(:,k);

        currentDiameters_av = (currentDiameters(1:end-1)+currentDiameters(2:end))/2; % This assumes that bins are only excluded due the instrument not having them at this stage
        currentVols = 4/3*pi*(currentDiameters_av/2).^3 * (1e-6)^3;

        if clipNegatives
            validSamples = currentData >= 0;

            maxExclusions = 3;
            m = 1;
            while ~validSamples(m) && m <= maxExclusions
                m = m+1;
            end

            currentData = currentData(m:end);
            currentData_raw = currentData_raw(m:end);
            currentLogBinSizes = currentLogBinSizes(m:end);
            currentDiameters = currentDiameters(m:end);
            currentDiameters_av = currentDiameters_av(m:end);
            currentVols = currentVols(m:end);

            currentData(currentData<0) = 0;
        end

        A_t = sum(currentData_raw);
        A_t_v = sum(currentData_raw .* currentVols);

        [~, mu_t, sigma_t, l_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(20), 'sig_LB', 0.2, 'sig_UB', 50, 'tubeCorrection', correctionVal);
        A_diff(k) = A_t;
        A_v_diff(k) = A_t_v;
        mu_diff(k) = mu_t;
        sigma_diff(k) = sigma_t;

    end

    % if clipNegatives
    %     A_fg_after(A_fg_after<0) = 0;
    %     %A_fg_after= A_fg_after(A_fg_after>0);
    % end
    
    
    fitLogNorm = true;
    if (fitLogNorm) % Fix should fit sum of lognorm and noise
        noiseMean = 0;
        noiseStd = 1117; %UG 1117. LG 3602
        
        reject = A_diff < -3*noiseStd; % Reject points that are probably errors
        
        A_diff_trunc_zeros = A_diff(A_diff>0);
        A_marg_ln_hat = lognfit(A_diff_trunc_zeros);
        A_diff_trunc = A_diff(~reject);
        objFun = @(x) -1*sum(sumLognormNormpdf(A_diff_trunc, x(1),x(2),0,noiseStd));

        x0 = A_marg_ln_hat;
        LB = [0,0];
        UB = [2*A_marg_ln_hat(1), 2*A_marg_ln_hat(2)];
        
        options_fmincon = optimoptions('fmincon','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, 'Display', 'off');
        bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
        A_marg_ln_hat = bestVal;
    else
        [A_marg_n_mu, A_marg_n_sig] = normfit(A_diff);
    end
    
    mu_weights = normcdf(A_diff, noiseMean+3*noiseStd,noiseStd);
    [mu_marg_mu, mu_marg_sig] = normfit(mu_diff,[],[],mu_weights);
    sig_marg_hat = gamfit(sigma_diff,[],[],mu_weights);
    
    mu_plot = linspace(-3,3,200);
    sig_plot = linspace(0,5,200);
    A_plot = linspace(-4e3,1e5,200);

    figure;
    subplot(3,2,1);
    scatter(exp(mu_diff),zeros(size(mu_diff)),'rx','LineWidth',2);
    xlabel('mean particle diameter (\mu m)');
    ylabel('density');
    title(['mu marginal: \mu = ', num2str(exp(mu_marg_mu)), ', \sigma = ', num2str(mu_marg_sig)]);
    hold on;
    plot(exp(mu_plot), normpdf(mu_plot,mu_marg_mu, mu_marg_sig), 'r');
    xlim([min(exp(mu_plot)), max(exp(mu_plot))]);

    subplot(3,2,3);
    scatter(sigma_diff,zeros(size(sigma_diff)),'bx','LineWidth',2);
    xlabel('sigma (log units)');
    ylabel('density');
    title(['sigma marginal: \alpha = ', num2str(sig_marg_hat(1)), ', \beta = ', num2str(sig_marg_hat(2))]);
    hold on;
    plot(sig_plot, gampdf(sig_plot,sig_marg_hat(1), sig_marg_hat(2)), 'b');
    xlim([min((sig_plot)), max((sig_plot))]);

    subplot(3,2,5);
    scatter(A_diff,zeros(size(A_diff)),'kx','LineWidth',2);
    xlabel('# particles/m^3/s');
    hold on;
    if fitLogNorm
        plot(A_plot, lognpdf(A_plot,A_marg_ln_hat(1),A_marg_ln_hat(2)), 'k');
        title(['#particles marginal: \mu = ', num2str(exp(A_marg_ln_hat(1))), ', \sigma = ', num2str(A_marg_ln_hat(2))]);
    else
        plot(A_plot, normpdf(A_plot,A_marg_n_mu,A_marg_n_sig), 'k');
        title(['#particles marginal: \mu = ', num2str(A_marg_n_mu), ', \sigma = ', num2str(A_marg_n_sig)]);
    end
    xlim([min(A_plot), max(A_plot)]);

    subplot(3,2,2);
    scatter(A_diff, exp(mu_diff),'gx','LineWidth',2);
    xlim([min(A_plot), max(A_plot)]);
    ylim([min(exp(mu_plot)), max(exp(mu_plot))]);
    xlabel('# particles/m^3/s');
    ylabel('particle diameter (\mu m)');
    title('Joint distribution \mu vs #');

    subplot(3,2,4);
    scatter(exp(mu_diff),sigma_diff,'gx','LineWidth',2);
    xlim([min(exp(mu_plot)), max(exp(mu_plot))]);
    ylim([min(sig_plot), max(sig_plot)]);
    xlabel('particle diameter (\mu m)');
    ylabel('sigma');
    title('Joint distribution \sigma vs \mu');

    subplot(3,2,6);
    scatter(A_diff,sigma_diff,'gx','LineWidth',2);
    xlim([min(A_plot), max(A_plot)]);
    ylim([min(sig_plot), max(sig_plot)]);
    xlabel('# particles/m^3/s');
    ylabel('sigma');
    title('Joint distribution \sigma vs #');
end

%% fit FG after the event
% Integrate to get in the same window
fgwindowSize = 50;
buffer = 20;
windowStart = -1; %Should really be 0 for all cases

sampleValid = time >= windowStart & time < (windowStart + fgwindowSize + buffer);
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

% if clipNegatives
%     sample_int_fg_after(sample_int_fg_after<0) = 0;
% end

figure;
for k=1:size(sample_int_fg_after,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    currentData = sample_int_fg_after(1:end,k);
    currentData_raw = currentData;
    currentDiameters = diameters_2(:,k);
    
    currentDiameters_av = (currentDiameters(1:end-1)+currentDiameters(2:end))/2; % This assumes that bins are only excluded due the instrument not having them at this stage
    currentVols = 4/3*pi*(currentDiameters_av/2).^3 * (1e-6)^3;
        
    if clipNegatives
        validSamples = currentData >= 0;
        
        m = 1;
        while ~validSamples(m) && m <= 2
            m = m+1;
        end
        
        currentData = currentData(m:end);
        currentData_raw = currentData_raw(m:end);
        currentLogBinSizes = currentLogBinSizes(m:end);
        currentDiameters = currentDiameters(m:end);
        currentDiameters_av = currentDiameters_av(m:end);
        currentVols = currentVols(m:end);
        
        currentData(currentData<0) = 0;
    end
    
    A_t = sum(currentData_raw);
    A_t_v = sum(currentData_raw .* currentVols);
    
    [~, mu_t, sigma_t, l_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(50), 'sig_LB', 0.2, 'sig_UB', 50);
    A_fg_after(k) = A_t;
    mu_fg_after(k) = mu_t;
    sigma_fg_after(k) = sigma_t;
    
    if k==4
        a = 1;
    end
    
end

%% fit FG before the event
buffer = 20;
windowEnd = windowStart; %Should really be 0 for all cases

sampleValid = time >= (windowStart- fgwindowSize - buffer) & time < (windowEnd);
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

figure;
for k=1:size(sample_int_fg_before,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    currentData = sample_int_fg_before(1:end,k);
    currentData_raw = currentData;
    currentDiameters = diameters_2(:,k);
    
    currentDiameters_av = (currentDiameters(1:end-1)+currentDiameters(2:end))/2; % This assumes that bins are only excluded due the instrument not having them at this stage
    currentVols = 4/3*pi*(currentDiameters_av/2).^3 * (1e-6)^3;
        
    if clipNegatives
        validSamples = currentData >= 0;
        
        m = 1;
        while ~validSamples(m) && m <= 2
            m = m+1;
        end
        
        currentData = currentData(m:end);
        currentData_raw = currentData_raw(m:end);
        currentLogBinSizes = currentLogBinSizes(m:end);
        currentDiameters = currentDiameters(m:end);
        currentDiameters_av = currentDiameters_av(m:end);
        currentVols = currentVols(m:end);
        
        currentData(currentData<0) = 0;
    end
    
    A_t = sum(currentData_raw);
    A_t_v = sum(currentData_raw .* currentVols);
    
    [~, mu_t, sigma_t, l_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(50), 'sig_LB', 0.2, 'sig_UB', 50);
    A_fg_before(k) = A_t;
    mu_fg_before(k) = mu_t;
    sigma_fg_before(k) = sigma_t;
    
    if k==7
        a = 1;
    end
    
end

if clipNegatives
    A_fg_before(A_fg_before<0) = 0;
    %A_fg_before = A_fg_before(A_fg_before>0);
end

%% fit BG after the event
% Integrate to get in the same window
windowSize = 80;
buffer = 20;
windowStart = -1; %Should really be 0 for all cases

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
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_bg_after(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(0.2), 'sig_LB', 0.1, 'sig_UB', 50);
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
    [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_bg_before(1:end,k)./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(0.2), 'sig_LB', 0.1, 'sig_UB', 50);
    A_bg_before(k) = A_t;
    mu_bg_before(k) = mu_t;
    sigma_bg_before(k) = sigma_t;
    
end

%% Now try to fit a bimodal distribution to the 'after' data
% Integrate to get in the same window
windowSize = 50;
buffer = 20;
windowStart = 0; %Should really be 0 for all cases

sampleValid = time >= windowStart & time < (windowStart + windowSize + buffer);
sample_preInt = raw(:,sampleValid,:);
sample_int_raw_after_all = zeros(size(raw,1),size(raw,3));

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
    sample_int_raw_after_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_raw_after_all(:,1)));

sample_int_raw_after = zeros(nValid,size(sample_int_raw_after_all,2));
diameters_2 = zeros(size(sample_int_raw_after,1)+1, size(sample_int_raw_after,2));

for k = 1:size(sample_int_bg_after_all,2)
    temp = sample_int_raw_after_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_raw_after(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

figure;
relative_mu_bounds = [-0.1,0.1];
relative_sig_bounds = [0.8, 1.2];
for k=1:size(sample_int_raw_after,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    
    mu_bg_est = mu_bg_after(k);
    sigma_bg_est = sigma_bg_after(k);
    %mu_fg_est = mu_fg_after(k);
    %sigma_fg_est = sigma_fg_after(k);
    mu_fg_est = mu_diff(k);
    sigma_fg_est = sigma_diff(k);
    
    A_1 = A_bg_after(k)/(normcdf(log(diameters_2(end,k)),mu_bg_est, sigma_bg_est) - normcdf(log(diameters_2(1,k)),mu_bg_est, sigma_bg_est));
    A_2 = A_fg_after(k)/(normcdf(log(diameters_2(end,k)),mu_fg_est, sigma_fg_est) - normcdf(log(diameters_2(1,k)),mu_fg_est, sigma_fg_est));
    w_est = abs(A_2) / abs(A_1);
    
    [A_t, mu_t, sigma_t, l_t, w_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_raw_after(1:end,k)./currentLogBinSizes)','fitType','counts', 'bimodal', true, 'mu_LB', mu_fg_est+relative_mu_bounds(1), 'mu_UB', min(mu_fg_est+relative_mu_bounds(2),4), 'sig_LB', sigma_fg_est*relative_sig_bounds(1), 'sig_UB', sigma_fg_est*relative_sig_bounds(2), 'bg_mu', mu_bg_est, 'bg_sig', sigma_bg_est, 'w_UB', w_est*100);
    A_after(k) = A_t;
    mu_after(k) = mu_t;
    sigma_after(k) = sigma_t;
    w_after(k) = w_t;
    
    if mu_t > 4
        a = 1;
    end
    
end

%% Fit bimodal to before data
buffer = 20;
windowEnd = windowStart; %Should really be 0 for all cases

sampleValid = time >= (windowStart- windowSize - buffer) & time < (windowEnd);
sample_preInt = fg(:,sampleValid,:);
sample_int_raw_before_all = zeros(size(fg,1),size(fg,3));

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
    sample_int_raw_before_all(:,k) = cumdd;
    
end

% Diameter array per each data set
nValid = nnz(~isnan(sample_int_raw_before_all(:,1)));

sample_int_raw_before = zeros(nValid,size(sample_int_raw_before_all,2));
diameters_2 = zeros(size(sample_int_raw_before,1)+1, size(sample_int_raw_before,2));

for k = 1:size(sample_int_bg_before_all,2)
    temp = sample_int_raw_before_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_raw_before(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end
for k = 1:size(sample_int_bg_before_all,2)
    temp = sample_int_raw_before_all(:,k);
    valid = ~isnan(temp);
    
    sample_int_raw_before(:,k) = temp(valid);
    diameters_2(:,k) = diameters([valid; true]);
end

figure;
relative_mu_bounds = [-0.1,0.1];
relative_sig_bounds = [0.95, 1.05];
for k=1:size(sample_int_raw_before,2)
    currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
    
    mu_bg_est = mu_bg_before(k);
    sigma_bg_est = sigma_bg_before(k);
    mu_fg_est = mu_after(k);
    sigma_fg_est = sigma_after(k);
    
    A_1 = A_bg_before(k)/(normcdf(log(diameters_2(end,k)),mu_bg_est, sigma_bg_est) - normcdf(log(diameters_2(1,k)),mu_bg_est, sigma_bg_est));
    A_2 = A_fg_before(k)/(normcdf(log(diameters_2(end,k)),mu_fg_est, sigma_fg_est) - normcdf(log(diameters_2(1,k)),mu_fg_est, sigma_fg_est));
    w_est = abs(A_2) / abs(A_1);
    
    [A_t, mu_t, sigma_t, l_t, w_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_raw_before(1:end,k)./currentLogBinSizes)','fitType','counts', 'bimodal', true, 'mu_LB', mu_fg_est+relative_mu_bounds(1), 'mu_UB', min(mu_fg_est+relative_mu_bounds(2),4), 'sig_LB', sigma_fg_est*relative_sig_bounds(1), 'sig_UB', sigma_fg_est*relative_sig_bounds(2), 'bg_mu', mu_bg_est, 'bg_sig', sigma_bg_est, 'w_UB', w_est*4);
    A_before(k) = A_t;
    mu_before(k) = mu_t;
    sigma_before(k) = sigma_t;
    w_before(k) = w_t;
    
    if mu_t > 4
        a = 1;
    end
    
end

figure;
scatter(w_before, w_after);
xlabel('w');
ylabel('w');
axis equal;

test = (w_after./w_before);

a = 1;

%% Plot
% subplot(3,2,1)
% plot(A_fg_before,A_fg_after,'x');
% hold on;
% plot(A_fg_before,A_fg_before);
% title('Amount of aerosol');
% 
% subplot(3,2,3)[A_t, mu_t, sigma_t, l_t, w_t] = fitAerosolDist(diameters_2(:,k).', (sample_int_raw_after(1:end,k)./currentLogBinSizes)','fitType','counts', 'bimodal', true, 'mu_LB', mu_fg_est-0.5, 'mu_UB', mu_fg_est+0.5, 'sig_LB', 0.9*sigma_fg_est, 'sig_UB', 1.1*sigma_fg_est, 'bg_mu', mu_bg_est, 'bg_sig', sigma_bg_est);
% plot(mu_fg_before,mu_fg_after,'x');
% hold on;
% plot(mu_fg_before,mu_fg_before);
% title('mu');
% 
% subplot(3,2,5)
% plot(sigma_fg_before,sigma_fg_after,'x');
% hold on;
% plot(sigma_fg_before,sigma_fg_before);
% title('sigma');
% 
% subplot(3,2,2)
% plot(A_bg_before,A_bg_after,'x');
% hold on;
% plot(A_bg_before,A_bg_before);
% title('Amount of aerosol');
% 
% subplot(3,2,4)
% plot(mu_bg_before,mu_bg_after,'x');
% hold on;
% plot(mu_bg_before,mu_bg_before);
% title('mu');
% 
% subplot(3,2,6)
% plot(sigma_bg_before,sigma_bg_after,'x');
% hold on;
% plot(sigma_bg_before,sigma_bg_before);
% title('sigma');

%% Fits
% [A_fg_before_mu, A_fg_before_sig] = normfit(A_fg_before);
% [A_fg_after_mu, A_fg_after_sig] = normfit(A_fg_after);
% [A_bg_before_mu, A_bg_before_sig] = normfit(A_bg_before);
% [A_bg_after_mu, A_bg_after_sig] = normfit(A_bg_after);
% A_fg_before_phat = gamfit(A_fg_before);
% A_fg_after_phat = gamfit(A_fg_after);
% A_bg_before_phat = gamfit(A_bg_before);
% A_bg_after_phat = gamfit(A_bg_after);
% A_fg_before_phat = gamfit(A_fg_before);
% A_fg_after_phat = gamfit(A_fg_after);
[A_fg_before_mu, A_fg_before_sig] = normfit(A_fg_before);
[A_fg_after_mu, A_fg_after_sig] = normfit(A_fg_after);
A_bg_before_phat = lognfit(A_bg_before);
A_bg_after_phat = lognfit(A_bg_after);

[mu_diff_mu, mu_diff_sig] = normfit(mu_diff);
[A_diff_mu, A_diff_sig] = normfit(A_diff);


[mu_fg_before_mu, mu_fg_before_sig] = normfit(mu_fg_before);
[mu_fg_after_mu, mu_fg_after_sig] = normfit(mu_fg_after);
[mu_bg_before_mu, mu_bg_before_sig] = normfit(mu_bg_before);
[mu_bg_after_mu, mu_bg_after_sig] = normfit(mu_bg_after);

%mu_fg_diff = mu_fg_after - mu_fg_before;
%mu_bg_diff = mu_bg_after - mu_bg_before;
%[mu_fg_diff_mu, mu_fg_diff_sig] = normfit(mu_fg_diff);
%[mu_bf_diff_mu, mu_bg_diff_sig] = normfit(mu_fg_diff);

sigma_fg_before_phat = gamfit(sigma_fg_before);
sigma_fg_after_phat = gamfit(sigma_fg_after);
sigma_bg_before_phat = gamfit(sigma_bg_before);
sigma_bg_after_phat = gamfit(sigma_bg_after);

%sig_
%sigma_fg_diff_phat = gamfit(sigma_fg_after - sig_bg_before);
%sigma_bg_diff_phat = gamfit(sigma_bg_after - sig_bg_before);

plot_A_bg = linspace(0,2e8,500);
plot_A_fg = linspace(-1e5,3e6,500);
plot_mu = linspace(-5,1,500);
plot_sig = linspace(0,2,500);


%rejectionSampleLognorm(sample_int_fg_after, diameters_2, mu_fg_after_mu, mu_fg_after_sig, sigma_fg_after_phat(1), sigma_fg_after_phat(2), mu_fg_after, sigma_fg_after)


figure;
subplot(6,1,1);
plot(exp(plot_mu), normpdf(plot_mu,mu_fg_before_mu, mu_fg_before_sig),'b');
hold on;
scatter(exp(mu_fg_before),zeros(size(mu_fg_before)),'xb','LineWidth',2);
xlim([min(exp(plot_mu)),max(exp(plot_mu))]);
title('mu fg dist');
xlabel('diameter (\mum)');
plot(exp(plot_mu), normpdf(plot_mu,mu_fg_after_mu, mu_fg_after_sig),'r');
hold on;
scatter(exp(mu_fg_after),ones(size(mu_fg_after))*1.1*max(normpdf(plot_mu,mu_fg_after_mu, mu_fg_after_sig)),'xr','LineWidth',2);
xlim([min(exp(plot_mu)),max(exp(plot_mu))]);
legend('before','','after');

subplot(6,1,2);
plot(exp(plot_sig), gampdf(plot_sig,sigma_fg_before_phat(1), sigma_fg_before_phat(2)),'b');
hold on;
scatter(exp(sigma_fg_before),zeros(size(sigma_fg_before)),'bx','LineWidth',2);
xlim([min(exp(plot_sig)),max(exp(plot_sig))]);
title('sigma fg dist');
plot(exp(plot_sig), gampdf(plot_sig,sigma_fg_after_phat(1), sigma_fg_after_phat(2)),'r');
hold on;
scatter(exp(sigma_fg_after),zeros(size(sigma_fg_after)),'rx','LineWidth',2);
xlim([min(exp(plot_sig)),max(exp(plot_sig))]);
xlabel('diameter (\mum)');
legend('before','','after');

subplot(6,1,3);
plot((plot_A_fg), normpdf(plot_A_fg,A_fg_before_mu, A_fg_before_sig), 'b');
%plot((plot_A_fg), gampdf(plot_A_fg,A_fg_before_phat(1), A_fg_before_phat(2)), 'b');
%plot((plot_A_fg), lognpdf(plot_A_fg,A_fg_before_phat(1), A_fg_before_phat(2)));
hold on;
scatter(A_fg_before,zeros(size(A_fg_before)),'bx','LineWidth',2);
plot((plot_A_fg), normpdf(plot_A_fg,A_fg_after_mu, A_fg_after_sig),'r');
%plot((plot_A_fg), gampdf(plot_A_fg,A_fg_after_phat(1), A_fg_after_phat(2)),'r');
%plot((plot_A_fg), lognpdf(plot_A_fg,A_fg_after_phat(1), A_fg_after_phat(2)));
scatter(A_fg_after,zeros(size(A_fg_after)),'rx','LineWidth',2);
%xlim([min(plot_A_fg),max(plot_A_fg)]);
title('A fg dist');
legend('before','','after');

subplot(6,1,4);
plot(exp(plot_mu), normpdf(plot_mu,mu_bg_before_mu, mu_bg_before_sig),'b');
hold on;
scatter(exp(mu_bg_before),zeros(size(mu_bg_before)),'bx','LineWidth',2);
%xlim([min(exp(plot_mu)),max(exp(plot_mu))]);
plot(exp(plot_mu), normpdf(plot_mu,mu_bg_after_mu, mu_bg_after_sig),'r');
scatter(exp(mu_bg_after),zeros(size(mu_bg_after)),'rx','LineWidth',2);
xlim([min(exp(plot_mu)),max(exp(plot_mu))]);
title('mu bg dist');
legend('before','','after');

subplot(6,1,5);
plot(exp(plot_sig), gampdf(plot_sig,sigma_bg_before_phat(1), sigma_bg_before_phat(2)),'b');
hold on;
scatter(exp(sigma_bg_before),zeros(size(sigma_bg_before)),'bx','LineWidth',2);
%xlim([min(exp(plot_sig)),max(exp(plot_sig))]);
plot(exp(plot_sig), gampdf(plot_sig,sigma_bg_after_phat(1), sigma_bg_after_phat(2)),'r');
hold on;
scatter(exp(sigma_bg_after),zeros(size(sigma_bg_after)),'rx','LineWidth',2);
xlim([min(exp(plot_sig)),max(exp(plot_sig))]);
title('sigma bg dist');
legend('before','','after');

subplot(6,1,6);
plot((plot_A_bg), lognpdf(plot_A_bg,A_bg_before_phat(1), A_bg_before_phat(2)),'b');
hold on;
scatter(A_bg_before,zeros(size(A_bg_before)),'bx','LineWidth',2);
%xlim([min(plot_A_fg),max(plot_A_fg)]);
plot((plot_A_bg), lognpdf(plot_A_bg,A_bg_after_phat(1), A_bg_after_phat(2)),'r');
hold on;
scatter(A_bg_after,zeros(size(A_bg_after)),'rx','LineWidth',2);
%xlim([min(plot_A_fg),max(plot_A_fg)]);
title('A bg dist');
legend('before','','after');


figure;
subplot(3,1,1);
scatter(mu_fg_before, mu_fg_after,'kx','LineWidth',2);
xlabel('mu before');
ylabel('mu after');
axis equal;

subplot(3,1,2);
scatter(A_fg_before,A_fg_after,'kx','LineWidth',2);
xlabel('A before');
ylabel('A after');
axis equal;

subplot(3,1,3);
scatter(sigma_fg_before,sigma_fg_after,'kx','LineWidth',2);
xlabel('sigma before');
ylabel('sigma after');
axis equal;


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

fileID = fopen('lognormal_inf.stan');
model_code = textscan(fileID,'%s');
model_code = model_code{:};
fclose(fileID);

%control.stepsize = 2;
%conrtol.stepsize_jitter = 10;
control.adapt_delta = 0.7;
control.adapt_kappa = 0.9;
%control.max_tree


sm = StanModel('model_code',model_code, 'model_name', 'lognormal_inf','verbose',true, 'init', initVals, 'chains', 1, 'iter', 10000, 'control', control, 'file_overwrite', true);
%sm.compile();

% subsequent calls will skip recompilation
fit = sm.sampling('data',stan_data);

%fit = stan('file','lognormal_inf.stan','data',stan_data,'verbose',true, 'init', initVals, 'chains', 4, 'iter', 4000);
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
