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
addpath('.\Violinplot-Matlab-master');

folder = 'C:\Users\george\OneDrive - The University of Nottingham\SAVE\csv0909_0411';

loadFolder = true;
if loadFolder
    fileList_raw = dir(folder);
    
    filterFun = @(x) regexpi(x, '.*[.]csv');
    temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
    fileList_raw = fileList_raw(~cellfun(@isempty,temp));
    
    filterFun = @(x) regexpi(x, '.*test.*[.]csv'); %Should exclude TEST files
    temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
    fileList_raw = fileList_raw(cellfun(@isempty,temp));
    
    filterFun = @(x) regexpi(x, '.*roomdoors.*[.]csv'); %Exclude room doors files
    temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
    fileList_raw = fileList_raw(cellfun(@isempty,temp));
    
    fileList = {fileList_raw.name};
    fileList = fileList';
    
    filterFun = @(x) regexp(x, '.*nullreference[.]csv');
    temp = cellfun(filterFun, fileList, 'UniformOutput', false); 
    fileList = fileList(cellfun(@isempty,temp));
else
    fileList = {'UpperGI_throatspraygiven.csv'};

end

upperGI_nullreffile = 'UpperGI_nullreference.csv';
lowerGI_nullreffile = 'LowerGI_nullreference.csv';

fileList = {upperGI_nullreffile, lowerGI_nullreffile, fileList{1:end}}';

nFiles = size(fileList, 1);
saveFigs = true;


%% Table for soring vars
varNames = {'label', 'Nevents', 'Npatients', 'mean_n', 'std_n', 'median_n', 'lq_n', 'uq_n',  'mean_v', 'std_v', 'median_v', 'lq_v', 'uq_v', 'mu_mean_n', 'mu_std_n', 'sig_a_n', 'sig_b_n', 'mu_mean_v', 'mu_std_v', 'sig_a_v', 'sig_b_v'};

for k=1:size(varNames,2)
    if k == 1
        varTypes{k} = 'string';
    else
        varTypes{k} = 'double';
    end
end

resultsTable = table('Size',[nFiles, size(varNames,2)],'VariableTypes', varTypes, 'VariableNames', varNames);

for fileIdx = 1:nFiles
    file = fileList{fileIdx};
    label = file(1:end-4);
    
    isLowerGI = ~isempty(regexp(file, 'LowerGI_*.'));
    isUpperGI = ~isempty(regexp(file, 'UpperGI_*.'));
    
    if fileIdx == 1
        if ~strcmpi(file, upperGI_nullreffile)
            error('Need null reference');
        end
        isUGI_nullref = true;
    else
        isUGI_nullref = false;
    end
    
    if fileIdx == 2
        if ~strcmpi(file, lowerGI_nullreffile)
            error('Need null reference');
        end
        isLGI_nullref = true;
    else
        isLGI_nullref = false;
    end
    
    filepath = fullfile(folder,file);

    % Load data
    opts = detectImportOptions(filepath);

    isFullData = true; %Not full data is for an earlier prototype version
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

        [tempTimes, tempIdxes] = unique(T(:,4));
        eventTimes = table2array(tempTimes);
        indices = 1:size(eventTimes,1);
        indices = indices';
        diameters_full = T.ParticleBin_um_;
        diameters = diameters_full(table2array(T(:,4)) == eventTimes(1));
        %diameters = table2array(diameters);
        dataStartCol = 25; %% FIX should determine this from variable names
        
        patientNos = T.StudyNumber(sort(tempIdxes, 'ascend'));
        
        otherVars = T(sort(tempIdxes, 'ascend'),8:(dataStartCol-2));
        
        % Sedation variable creations
        sedation_raw = otherVars.SedationDetails;
        [sedation, sedationCats] = processVars(otherVars.SedationDetails, {'midazolam','propofol', 'entonox', 'none'}, {'.*mid.*'}, {'.*prop.*'}, {'.*ento.*'},{'.*none.*','.*N/A.*',''});
                
        % Anal tone
        analTone_raw = otherVars.AnalTone;
        [analTone, analToneCats] = processVars(otherVars.AnalTone, {'low','medium', 'high', 'N/A'}, {'.*low.*', '.*loose.*', '.*poor.*'}, {'.*normal.*', '.*medium.*', '.*middle.*'}, {'.*good.*', '.*tight.*', '.*high.*'},{'.*none.*','.*N/A.*',''});
                
        %Co2 vs water
        co2water_raw = otherVars.UseOfCO2OrWater;
        [useOfCO2orWater, useOfCO2orWaterCats] = processVars(otherVars.UseOfCO2OrWater, {'CO2','Water', 'N/A'}, {'.*co2.*', '.*yes.*'}, {'.*water.*', '.*h2o.*'}, {'.*none.*','.*N/A.*',''});
             
        % Smoker
        smokerTemp_raw = otherVars.Smoker;
        [isSmoker, isSmokerCats] = processVars(otherVars.Smoker, {'yes','no', 'unknown'}, {'.*yes.*', '.*[1-9]?[0-9]*.*'}, {'.*no.*', '.*[0]+.*'}, {'.*none.*','.*N/A.*',''});
        
        % Mask
        usesPatientMask_raw = ismember(patientNos,[59,60,61, 64,65]); %FIX should pull this
        usesPatientMask = categorical(usesPatientMask_raw+1,1:2,{'No Mask', 'Mask'});
        % Room type
        
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


    nEvents = size(indices,1);
    nPatients = size(unique(patientNos),1);

    maxDiameter = 25; % Should load this from datasheet?
    diameters = [diameters; maxDiameter];

    time = startTime:timeStep:endTime;
    nTimes = (endTime - startTime)/timeStep + 1;

    validRows = table2array(T(:,4)) == eventTimes(1); % FIX should loop for other indices and check
    tempData = table2array(T(validRows,dataStartCol:dataStartCol+nTimes-1));
    bg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
    fg = zeros(size(tempData,1), size(tempData,2), size(indices,1));
    raw = zeros(size(tempData,1), size(tempData,2), size(indices,1));
    avSampleTimes = zeros(size(indices,1),1);

    %% Plot raw data
    if fileIdx == 1
        rawDataFig = figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure(rawDataFig);
    end
    clf;
    for currentIdxIdx = 1:size(indices,1)
        
        currentIdx = indices(currentIdxIdx);
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

        tempValid = ~isnan(data);
        tempValid = nansum(tempValid,1) > 0;
        avSampleTime = median(diff(time(tempValid)));
        avSampleTimes(currentIdx) = avSampleTime;
        
        correctForAT = true; %AeroTrak time is the start of the sample
        if correctForAT
            data = circshift(data, avSampleTime);
            data(1:avSampleTime) = NaN;
            data_v = circshift(data_v, avSampleTime);
            data_v(1:avSampleTime) = NaN;
        end
        
        % Densities so that a probability density approach can be used
        data_v_density = data_v ./ repmat(currentLogBinSizes,1, size(data,2)); %Try using log binsizes
        data_density = data ./ repmat(currentLogBinSizes,1, size(data,2));
        
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

    if saveFigs
        saveFigName = [label, '_superpos'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end

    clipNegatives = true; %Negative counts values set to zero


    %% Raw difference
    calcRawRatio = true;

    if calcRawRatio
        rawratiowindowSize_after = 100;
        buffer = 20;
        windowStart = 0; %Should really be 0 for all cases

        sampleValid = time > windowStart & time <= (windowStart + rawratiowindowSize_after + buffer);
        sample_preInt = raw(:,sampleValid,:);
        sample_int_raw_after_all = zeros(size(raw,1),size(raw,3));

        for k = 1:size(sample_preInt,3)
            cumdt = 0;
            cumdd = zeros(size(sample_preInt,1),1);
            realWindowSize = 0;
                   
            for kk=1:size(sample_preInt,2)

                if (kk <= rawratiowindowSize_after)
                    cumdt = cumdt + 1;
                end
                
                if cumdt > avSampleTimes(k)
                    cumdd = cumdd - lastSample;
                    realWindowSize = realWindowSize - avSampleTimes(k); %In case there aren't enough afterwards
                    break;
                end

                if ~all(isnan(sample_preInt(:,kk,k)))
                    lastSample = sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);
                    cumdd = cumdd + lastSample;

                    realWindowSize = realWindowSize+cumdt;
                    cumdt = 0;

                    if kk >= rawratiowindowSize_after
                        break;
                    end
                end

            end

            cumdd = cumdd/realWindowSize;
            sample_int_raw_after_all(:,k) = cumdd;
        end

        %Integrate before
        windowEnd = windowStart;
        rawratiowindowSize_before = 100;
        windowStart = windowEnd - rawratiowindowSize_before;

        sampleValid = time > windowStart & time <= (windowStart + rawratiowindowSize_before + buffer);
        sample_preInt = raw(:,sampleValid,:);
        sample_int_raw_before_all = zeros(size(raw,1),size(raw,3));

        for k = 1:size(sample_preInt,3)
            cumdt = 0;
            cumdd = zeros(size(sample_preInt,1),1);
            realWindowSize = 0;
            for kk=1:size(sample_preInt,2)

                if (kk <= rawratiowindowSize_before)
                    cumdt = cumdt + 1;
                end

                if ~all(isnan(sample_preInt(:,kk,k)))
                    if cumdt > avSampleTimes(k)
                        cumdt = avSampleTimes(k); % In case there aren't enough samples before
                    end
                    
                    cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);

                    realWindowSize = realWindowSize+cumdt;
                    cumdt = 0;

                    if kk >= rawratiowindowSize_before
                        break;
                    end
                end

            end

            cumdd = cumdd/realWindowSize;
            sample_int_raw_before_all(:,k) = cumdd;

        end

        sample_int_raw_ratio_all = sample_int_raw_after_all./sample_int_raw_before_all;

        % Diameter array per each data set
        nValid = nnz(~isnan(sample_int_raw_ratio_all(:,end))); %FIX Should be more robust than this

        sample_int_raw_ratio = zeros(nValid,size(sample_int_raw_ratio_all,2));
        diameters_2 = zeros(size(sample_int_raw_ratio,1)+1, size(sample_int_raw_ratio,2));

        for k = 1:size(sample_int_raw_ratio_all,2)
            temp = sample_int_raw_ratio_all(:,k);
            valid = ~isnan(temp);

            sample_int_raw_ratio(:,k) = temp(valid);
            diameters_2(:,k) = diameters([valid; true]);
        end


        if fileIdx == 1
            curveFitFig = figure;
        else
            figure(curveFitFig);
        end
      
        A_ratio = [];
        A_v_ratio = [];
        mu_ratio = [];
        sigma_ratio = [];
        mu_v_ratio = [];
        sigma_v_ratio = [];
        
        for k=1:size(sample_int_raw_ratio,2)
            %figure;
            currentLogBinSizes = log(diameters_2(2:end,k)) - log(diameters_2(1:end-1,k));
            currentData = sample_int_raw_ratio(1:end,k);
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

            [A_est, mu_t, sigma_t, l_t, ~, A_v_est, mu_v_t, sigma_v_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(20), 'sig_LB', 0.2, 'sig_UB', 50, 'tubeCorrection', correctionVal, 'fullDiameters', diameters_2(:,k));
           
            % Can't use A_est as it artificially inflates larger numbers
            %if A_est == 0
                A_ratio(k) = A_t;
                A_v_ratio(k) = A_t_v;
            %else
            %    A_diff(k) = A_est;
            %    A_v_diff(k) = A_v_est;
            %end
     
            mu_ratio(k) = mu_t;
            sigma_diff(k) = sigma_t;
            mu_v_ratio(k) = mu_v_t;
            sigma_v_diff(k) = sigma_v_t;
            
            directFit_v = false;
            if directFit_v
                [A_v_est2, mu_v_t2, sigma_v_t2] = fitAerosolDist(currentDiameters.', (currentData.*currentVols./currentLogBinSizes)','fitType','volume', 'mu_LB', log(0.1), 'mu_UB', log(20), 'sig_LB', 0.2, 'sig_UB', 50, 'tubeCorrection', correctionVal, 'fullDiameters', diameters_2(:,k));
                mu_v_ratio(k) = mu_v_t2;
                sigma_v_diff(k) = sigma_v_t2;
            end

        end

        % if clipNegatives
        %     A_fg_after(A_fg_after<0) = 0;
        %     %A_fg_after= A_fg_after(A_fg_after>0);
        % end

        
        if (isUGI_nullref)
            [A_marg_n_mu, A_marg_n_sig] = normfit(A_ratio);
            noiseMean_UG = 0;
            noiseStd_UG = A_marg_n_sig;
            
            [A_v_marg_n_mu, A_v_marg_n_sig] = normfit(A_v_ratio);
            noiseMean_v_UG = 0;
            noiseStd_v_UG = A_v_marg_n_sig;
            
            noiseMean = noiseMean_UG;
            noiseStd = noiseStd_UG;
            noiseMean_v = noiseMean_v_UG;
            noiseStd_v = noiseStd_v_UG;
        elseif isLGI_nullref
            [A_marg_n_mu, A_marg_n_sig] = normfit(A_ratio);
            noiseMean_LG = 0;
            noiseStd_LG = A_marg_n_sig;
            
            [A_v_marg_n_mu, A_v_marg_n_sig] = normfit(A_v_ratio);
            noiseMean_v_LG = 0;
            noiseStd_v_LG = A_v_marg_n_sig;
            
            
            noiseMean = noiseMean_LG;
            noiseStd = noiseStd_LG;
            noiseMean_v = noiseMean_v_LG;
            noiseStd_v = noiseStd_v_LG;
        else
            fitLogNorm = true;
            if (fitLogNorm) % Fix should fit sum of lognorm and noise
                
                if isUpperGI
                    noiseMean = noiseMean_UG;
                    noiseStd = noiseStd_UG;
                    noiseMean_v = noiseMean_v_UG;
                    noiseStd_v = noiseStd_v_UG;
                elseif isLowerGI
                    noiseMean = noiseMean_LG;
                    noiseStd = noiseStd_LG;
                    noiseMean_v = noiseMean_v_LG;
                    noiseStd_v = noiseStd_v_LG;
                else
                    error('Can''t determine noise reference!');
                end

                % Fir number
                reject = A_ratio < -3*noiseStd; % Reject points that are probably errors
                
                if nnz(~reject) == 0
                    reject = false(size(reject));
                end

                A_diff_trunc_zeros = A_ratio(A_ratio>0);
                A_marg_ln_hat = lognfit(A_diff_trunc_zeros);
                A_diff_trunc = A_ratio(~reject);
                objFun = @(x) -1*sum(sumLognormNormpdf(A_diff_trunc, x(1),x(2),0,noiseStd));
                
                if isnan(A_marg_ln_hat(1))
                    A_marg_ln_hat(1) = 7;
                end
                  
                if isnan(A_marg_ln_hat(2))
                    A_marg_ln_hat(2) = 7;
                end

                x0 = A_marg_ln_hat;
                LB = [0,0];
                UB = [2*A_marg_ln_hat(1), 2*A_marg_ln_hat(2)];

                options_fmincon = optimoptions('fmincon','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, 'Display', 'off');
                bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
                A_marg_ln_hat = bestVal;
                
                % Fit vol
                reject = A_v_ratio < -3*noiseStd_v; % Reject points that are probably errors
                
                if nnz(~reject) == 0
                    reject = false(size(reject));
                end

                A_v_diff_trunc_zeros = A_v_ratio(A_v_ratio>0);
                A_v_marg_ln_hat = lognfit(A_v_diff_trunc_zeros);
                A_v_diff_trunc = A_v_ratio(~reject);
                objFun = @(x) -1*sum(sumLognormNormpdf(A_v_diff_trunc, x(1),x(2),0,noiseStd_v));
                
                if isnan(A_v_marg_ln_hat(1))
                    A_v_marg_ln_hat(1) = 7;
                end
                  
                if isnan(A_v_marg_ln_hat(2))
                    A_v_marg_ln_hat(2) = 7;
                end

                x0 = A_v_marg_ln_hat;
                LB = [0,0];
                UB = [2*A_v_marg_ln_hat(1), 2*A_v_marg_ln_hat(2)];

                options_fmincon = optimoptions('fmincon','MaxFunctionEvaluations',1e4, 'StepTolerance',1e-12, 'OptimalityTolerance',1e-12, 'Display', 'off');
                bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
                A_v_marg_ln_hat = bestVal;
            else
                [A_marg_n_mu, A_marg_n_sig] = normfit(A_ratio);
            end
        end

        if (~isUGI_nullref && ~isLGI_nullref)
            mu_weights = normcdf(A_ratio, noiseMean+3*noiseStd,noiseStd);
            mu_weights = mu_weights/sum(mu_weights)*10;
            [mu_marg_mu, mu_marg_sig] = normfit(mu_ratio,[],[],mu_weights);
            sig_marg_hat = gamfit(sigma_diff,[],[],mu_weights);
            
            mu_v_weights = normcdf(A_v_ratio, noiseMean_v+3*noiseStd_v,noiseStd_v);
            mu_v_weights = mu_v_weights/sum(mu_v_weights)*10;
            [mu_v_marg_mu, mu_v_marg_sig] = normfit(mu_v_ratio,[],[],mu_v_weights*100);
            sig_v_marg_hat = gamfit(sigma_v_diff,[],[],mu_v_weights*100);
        else
            [mu_marg_mu, mu_marg_sig] = normfit(mu_ratio);
            sig_marg_hat = gamfit(sigma_diff);
            
            [mu_v_marg_mu, mu_v_marg_sig] = normfit(mu_v_ratio);
            sig_v_marg_hat = gamfit(sigma_v_diff);
        end

        mu_plot = linspace(-3,3,200);
        sig_plot = linspace(0,5,200);
        A_plot = linspace(-4e3,1e5,200);

        if (fileIdx == 1)
            nStatsFig = figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure(nStatsFig);
        end
        clf;
        subplot(3,2,1);
        scatter(exp(mu_ratio),zeros(size(mu_ratio)),'rx','LineWidth',2);
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
        scatter(A_ratio,zeros(size(A_ratio)),'kx','LineWidth',2);
        xlabel('# particles/m^3/s');
        hold on;
        
        if (isUGI_nullref || isLGI_nullref)
            plot(A_plot, normpdf(A_plot,A_marg_n_mu,A_marg_n_sig), 'k');
            title(['#particles marginal: \mu = ', num2str(A_marg_n_mu), ', \sigma = ', num2str(A_marg_n_sig)]);
        else
            if fitLogNorm
                plot(A_plot, lognpdf(A_plot,A_marg_ln_hat(1),A_marg_ln_hat(2)), 'k');
                title(['#particles marginal: \mu = ', num2str(exp(A_marg_ln_hat(1))), ', \sigma = ', num2str(A_marg_ln_hat(2))]);
            else
                plot(A_plot, normpdf(A_plot,A_marg_n_mu,A_marg_n_sig), 'k');
                title(['#particles marginal: \mu = ', num2str(A_marg_n_mu), ', \sigma = ', num2str(A_marg_n_sig)]);
            end
        end
        xlim([min(A_plot), max(A_plot)]);

        subplot(3,2,2);
        scatter(A_ratio, exp(mu_ratio),'gx','LineWidth',2);
        xlim([min(A_plot), max(A_plot)]);
        ylim([min(exp(mu_plot)), max(exp(mu_plot))]);
        xlabel('# particles/m^3/s');
        ylabel('particle diameter (\mu m)');
        title('Joint distribution \mu vs #');

        subplot(3,2,4);
        scatter(exp(mu_ratio),sigma_diff,'gx','LineWidth',2);
        xlim([min(exp(mu_plot)), max(exp(mu_plot))]);
        ylim([min(sig_plot), max(sig_plot)]);
        xlabel('particle diameter (\mu m)');
        ylabel('sigma');
        title('Joint distribution \sigma vs \mu');

        subplot(3,2,6);
        scatter(A_ratio,sigma_diff,'gx','LineWidth',2);
        xlim([min(A_plot), max(A_plot)]);
        ylim([min(sig_plot), max(sig_plot)]);
        xlabel('# particles/m^3/s');
        ylabel('sigma');
        title('Joint distribution \sigma vs #');
        
        if (saveFigs)
            saveas(gcf,fullfile(folder, [label,'_stats.fig']));
            saveas(gcf,fullfile(folder, [label,'_stats.png']));
        end
        
        %% Now plot volumes    
        A_v_plot = linspace(-2e-13,5e-12,200);

        if (fileIdx == 1)
            vStatsFig = figure('units','normalized','outerposition',[0 0 1 1]);
        else
            figure(vStatsFig);
        end
        clf;
        subplot(3,2,1);
        scatter(exp(mu_v_ratio),zeros(size(mu_v_ratio)),'rx','LineWidth',2);
        xlabel('mean particle diameter (\mu m)');
        ylabel('density');
        title(['mu marginal: \mu = ', num2str(exp(mu_marg_mu)), ', \sigma = ', num2str(mu_marg_sig)]);
        hold on;
        plot(exp(mu_plot), normpdf(mu_plot,mu_v_marg_mu, mu_v_marg_sig), 'r');
        xlim([min(exp(mu_plot)), max(exp(mu_plot))]);

        subplot(3,2,3);
        scatter(sigma_v_diff,zeros(size(sigma_v_diff)),'bx','LineWidth',2);
        xlabel('sigma (log units)');
        ylabel('density');
        title(['sigma marginal: \alpha = ', num2str(sig_v_marg_hat(1)), ', \beta = ', num2str(sig_v_marg_hat(2))]);
        hold on;
        plot(sig_plot, gampdf(sig_plot,sig_v_marg_hat(1), sig_v_marg_hat(2)), 'b');
        xlim([min((sig_plot)), max((sig_plot))]);

        subplot(3,2,5);
        scatter(A_v_ratio,zeros(size(A_v_ratio)),'kx','LineWidth',2);
        xlabel('vol particles/m^3/s');
        hold on;
        
        if (isUGI_nullref || isLGI_nullref)
            plot(A_v_plot, normpdf(A_v_plot,A_v_marg_n_mu,A_v_marg_n_sig), 'k');
            title(['vol particles marginal: \mu = ', num2str(A_v_marg_n_mu), ', \sigma = ', num2str(A_v_marg_n_sig)]);
        else
            if fitLogNorm
                plot(A_v_plot, lognpdf(A_v_plot,A_v_marg_ln_hat(1),A_v_marg_ln_hat(2)), 'k');
                title(['vol particles marginal: \mu = ', num2str(exp(A_v_marg_ln_hat(1))), ', \sigma = ', num2str(A_v_marg_ln_hat(2))]);
            else
                plot(A_v_plot, normpdf(A_v_plot,A_v_marg_n_mu,A_v_marg_n_sig), 'k');
                title(['vol particles marginal: \mu = ', num2str(A_v_marg_n_mu), ', \sigma = ', num2str(A_v_marg_n_sig)]);
            end
        end
        xlim([min(A_v_plot), max(A_v_plot)]);

        subplot(3,2,2);
        scatter(A_v_ratio, exp(mu_v_ratio),'gx','LineWidth',2);
        xlim([min(A_v_plot), max(A_v_plot)]);
        ylim([min(exp(mu_plot)), max(exp(mu_plot))]);
        xlabel('vol particles/m^3/s');
        ylabel('particle diameter (\mu m)');
        title('Joint distribution \mu vs #');

        subplot(3,2,4);
        scatter(exp(mu_v_ratio),sigma_v_diff,'gx','LineWidth',2);
        xlim([min(exp(mu_plot)), max(exp(mu_plot))]);
        ylim([min(sig_plot), max(sig_plot)]);
        xlabel('particle diameter (\mu m)');
        ylabel('sigma');
        title('Joint distribution \sigma vs \mu');

        subplot(3,2,6);
        scatter(A_v_ratio,sigma_v_diff,'gx','LineWidth',2);
        xlim([min(A_v_plot), max(A_v_plot)]);
        ylim([min(sig_plot), max(sig_plot)]);
        xlabel('vol particles/m^3/s');
        ylabel('sigma');
        title('Joint distribution \sigma vs #');
        
        if saveFigs
            saveas(gcf,fullfile(folder, [label,'_v_stats.fig']));
            saveas(gcf,fullfile(folder, [label,'_v_stats.png']));
        end
    end
    
    %% Now compile table
    resultsTable.label(fileIdx) = label;
    resultsTable.Nevents(fileIdx) = nEvents;
    resultsTable.Npatients(fileIdx) = nPatients;
    if (isUGI_nullref || isLGI_nullref)
        resultsTable.mean_n(fileIdx) = log(A_marg_n_mu);
        resultsTable.std_n(fileIdx) = log(A_marg_n_sig);
        resultsTable.mean_v(fileIdx) = log(A_v_marg_n_mu);
        resultsTable.std_v(fileIdx) = log(A_v_marg_n_sig);
    else
        if fitLogNorm
            resultsTable.mean_n(fileIdx) = A_marg_ln_hat(1);
            resultsTable.std_n(fileIdx) = A_marg_ln_hat(2);
            resultsTable.mean_v(fileIdx) = A_v_marg_ln_hat(1);
            resultsTable.std_v(fileIdx) = A_v_marg_ln_hat(2);
        else
            resultsTable.mean_n(fileIdx) = A_marg_n_mu;
            resultsTable.std_n(fileIdx) = A_marg_n_sigma;
            resultsTable.mean_v(fileIdx) = A_v_marg_n_mu;
            resultsTable.std_v(fileIdx) = A_v_marg_n_sigma;
        end
    end
    
    resultsTable.median_n(fileIdx) = median(A_ratio);
    resultsTable.lq_n(fileIdx) = quantile(A_ratio, 0.25);
    resultsTable.uq_n(fileIdx) = quantile(A_ratio, 0.75);
      
    resultsTable.median_v(fileIdx) = median(A_v_ratio);
    resultsTable.lq_v(fileIdx) = quantile(A_v_ratio, 0.25);
    resultsTable.uq_v(fileIdx) = quantile(A_v_ratio, 0.75);
    
    resultsTable.mu_mean_n(fileIdx) = mu_marg_mu;
    resultsTable.mu_std_n(fileIdx) = mu_marg_sig;  
    resultsTable.mu_mean_v(fileIdx) = mu_v_marg_mu;
    resultsTable.mu_std_v(fileIdx) = mu_v_marg_sig;
      
    resultsTable.sig_a_n(fileIdx) = sig_marg_hat(1);
    resultsTable.sig_b_n(fileIdx) = sig_marg_hat(2);
    resultsTable.sig_a_v(fileIdx) = sig_v_marg_hat(1);
    resultsTable.sig_b_v(fileIdx) = sig_v_marg_hat(2);
    
    resultsTable.n_raw(fileIdx) = {A_ratio};
    resultsTable.v_raw(fileIdx) = {A_v_ratio};
    resultsTable.mu_raw(fileIdx) = {mu_ratio};
    resultsTable.mu_v_raw(fileIdx) = {mu_v_ratio};
    
    resultsTable.age(fileIdx) = {otherVars.Age};
    resultsTable.bmi(fileIdx) = {otherVars.BMI};
    resultsTable.sedation(fileIdx) = {sedation};
    resultsTable.analTone(fileIdx) = {analTone};
    resultsTable.useOfCO2OrWater(fileIdx) = {useOfCO2orWater};
    resultsTable.isSmoker(fileIdx) = {isSmoker};
    resultsTable.usesPatientMask(fileIdx) = {usesPatientMask};
    resultsTable.patientNos(fileIdx) = {patientNos};
    
    data_box = [];
    data_v_box = [];
    mu_box = [];
    mu_v_box = [];
    cats_box = [];
    for dataSetIdx = 1:size(resultsTable,1)
        currentRaw_n = resultsTable.n_raw{dataSetIdx};
        currentRaw_n = currentRaw_n(:);
        currentRaw_v = resultsTable.v_raw{dataSetIdx};
        currentRaw_v = currentRaw_v(:);
        
        currentRaw_mu_n = resultsTable.mu_raw{dataSetIdx};
        currentRaw_mu_n = currentRaw_mu_n(:);
        currentRaw_mu_v = resultsTable.mu_v_raw{dataSetIdx};
        currentRaw_mu_v = currentRaw_mu_v(:);
        
        current_cats = repmat(resultsTable.label(dataSetIdx),size(currentRaw_n,1),1);
        current_cats = categorical(current_cats);
        
        data_box = [data_box; currentRaw_n];
        data_v_box = [data_v_box; currentRaw_v];
        
        temp = resultsTable.n_raw(dataSetIdx);
        if (max(temp{1}) > noiseMean + 3 *noiseStd) % Don't plot if not statistically significant
            mu_box = [mu_box; currentRaw_mu_n;];
            mu_v_box = [mu_v_box; currentRaw_mu_v;];
        else
            mu_box = [mu_box; zeros(size(currentRaw_mu_n));];
            mu_v_box = [mu_v_box; zeros(size(currentRaw_mu_v));];
        end
            
        cats_box = [cats_box; current_cats;];
    end
    
    if (fileIdx == 1)
        boxFig_n = figure('units','normalized','outerposition',[0 0 1 1]);
        violinFig_n = figure('units','normalized','outerposition',[0 0 1 1]);
        boxFig_v = figure('units','normalized','outerposition',[0 0 1 1]);
        violinFig_v = figure('units','normalized','outerposition',[0 0 1 1]);
        
        boxFig_mu_n = figure('units','normalized','outerposition',[0 0 1 1]);
        violinFig_mu_n = figure('units','normalized','outerposition',[0 0 1 1]);
        boxFig_mu_v = figure('units','normalized','outerposition',[0 0 1 1]);
        violinFig_mu_v = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    
    figure(boxFig_n);
    clf;
    boxplot(data_box,cats_box);
    tempCats_ = unique(cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    hold on;
    scatter(1:fileIdx,exp(resultsTable.mean_n(1:fileIdx)), 150,'kx', 'LineWidth', 3); 
    ylabel('Number of particles/m^3/s');
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(data_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht,'Rotation',60)
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['particle_number_box'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    idx1 = 1;
    idx2 = 2;
    temp_n = resultsTable.n_raw{idx1};
    dataSet1 = temp_n(:);
    temp_n = resultsTable.n_raw{idx2};
    dataSet2 = temp_n(:);
    %computeSignificance(dataSet1, dataSet2, noiseMean, noiseStd, resultsTable.mean_n(idx1), resultsTable.std_n(idx1), resultsTable.mean_n(idx2), resultsTable.std_n(idx2));

    
    figure(violinFig_n);
    clf;
    violinplot(data_box,cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    hold on;
    scatter(1:fileIdx,exp(resultsTable.mean_n(1:fileIdx)), 150,'kx', 'LineWidth', 3);
    ylabel('Number of particles/m^3/s');
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(data_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['particle_number_violin'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(boxFig_v);
    clf;
    boxplot(data_v_box,cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    hold on;
    scatter(1:fileIdx,exp(resultsTable.mean_v(1:fileIdx)), 150,'kx', 'LineWidth', 3);
    ylabel('Volume of particles/m^3/s');
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(data_v_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['particle_volume_box'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(violinFig_v);
    clf;
    violinplot(data_v_box,cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    hold on;
    ylim_hold = ylim();
    scatter(1:fileIdx,exp(resultsTable.mean_v(1:fileIdx)), 150,'kx', 'LineWidth', 3);
    ylim(ylim_hold);
    ylabel('Volume of particles/m^3/s');
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(data_v_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['particle_volume_violin'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(boxFig_mu_n);
    clf;
    boxplot(mu_box,cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    ylabel('Mean particle diameter (\mum)');
    y_ticks = yticklabels;
    for yidx = 1:size(y_ticks,1)
        temp = str2num(y_ticks{yidx});
        temp = exp(temp);
        temp = sprintf('%2.2f',temp);
        y_ticks{yidx} = temp;
    end
    yticklabels(y_ticks);
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(mu_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['diameter_number_box'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(violinFig_mu_n);
    clf;
    violinplot(mu_box,cats_box);
    set(gca,'TickLabelInterpreter', 'none');
    xtickangle(60);
    ylabel('Mean particle diameter (\mum)');
    y_ticks = yticklabels;
    for yidx = 1:size(y_ticks,1)
        temp = str2num(y_ticks{yidx});
        temp = exp(temp);
        temp = sprintf('%2.2f',temp);
        y_ticks{yidx} = temp;
    end
    yticklabels(y_ticks);
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(mu_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['diameter_number_violin'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(boxFig_mu_v);
    clf;
    boxplot(mu_v_box,cats_box);
    xtickangle(60);
    set(gca,'TickLabelInterpreter', 'none');
    ylabel('Mean particle diameter, volume fit(\mum)');
    y_ticks = yticklabels;
    for yidx = 1:size(y_ticks,1)
        temp = str2num(y_ticks{yidx});
        temp = exp(temp);
        temp = sprintf('%2.2f',temp);
        y_ticks{yidx} = temp;
    end
    yticklabels(y_ticks);
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(mu_v_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['diameter_volume_box'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(violinFig_mu_v);
    clf;
    violinplot(mu_v_box,cats_box);
    xtickangle(60);
    set(gca,'TickLabelInterpreter', 'none');
    ylabel('Mean particle diameter, volume fit(\mum)');
    y_ticks = yticklabels;
    for yidx = 1:size(y_ticks,1)
        temp = str2num(y_ticks{yidx});
        temp1 = exp(temp);
        temp2 = sprintf('%2.2f',temp1);
        y_ticks{yidx} = temp2;
    end
    yticklabels(y_ticks);
    for tempIdx = 1:fileIdx
        currentNevent = resultsTable.Nevents(tempIdx);
        currentNpatient = resultsTable.Npatients(tempIdx);
        ht = text(tempIdx,max(mu_v_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = ['diameter_volume_violin'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    if (fileIdx == 1)
        sedationFig_n = figure;
        ageFig_n = figure;
        bmiFig_n = figure;
        analToneFig_n = figure;
        co2vWaterFig_n = figure;
        smokerFig_n = figure;
        maskFig_n = figure;
       
        %repeat for volume
    end
    
    figure(sedationFig_n);
    clf;
    temp = resultsTable.sedation(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    boxplot(temp_n, temp);
    ylabel('Number of particles/m^3/s');
    xlabel('Sedation');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    tempCats = unique(temp);
    nCats = size(tempCats,1);
    for catIdx = 1:nCats
        currentNevent = nnz(temp == tempCats(catIdx));
        currentVals = temp_n(temp == tempCats(catIdx));
        temp_p_no = resultsTable.patientNos(fileIdx);
        temp_p_no = temp_p_no{1};
        currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
        ht = text(catIdx,max(currentVals)*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = [label, '_sedation'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(analToneFig_n);
    clf;
    temp = resultsTable.analTone(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    boxplot(temp_n, temp);
    ylabel('Number of particles/m^3/s');
    xlabel('Anal tone');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    tempCats = unique(temp);
    nCats = size(tempCats,1);
    for catIdx = 1:nCats
        currentNevent = nnz(temp == tempCats(catIdx));
        currentVals = temp_n(temp == tempCats(catIdx));
        temp_p_no = resultsTable.patientNos(fileIdx);
        temp_p_no = temp_p_no{1};
        currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
        ht = text(catIdx,max(currentVals)*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = [label, '_analtone'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(co2vWaterFig_n);
    clf;
    temp = resultsTable.useOfCO2OrWater(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    boxplot(temp_n, temp);
    ylabel('Number of particles/m^3/s');
    xlabel('Use of CO_2 vs. water');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    tempCats = unique(temp);
    nCats = size(tempCats,1);
    for catIdx = 1:nCats
        currentNevent = nnz(temp == tempCats(catIdx));
        currentVals = temp_n(temp == tempCats(catIdx));
        temp_p_no = resultsTable.patientNos(fileIdx);
        temp_p_no = temp_p_no{1};
        currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
        ht = text(catIdx,max(currentVals)*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = [label, '_CO2orWater'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(smokerFig_n);
    clf;
    temp = resultsTable.isSmoker(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    boxplot(temp_n, temp);
    ylabel('Number of particles/m^3/s');
    xlabel('Smoker');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    tempCats = unique(temp);
    nCats = size(tempCats,1);
    for catIdx = 1:nCats
        currentNevent = nnz(temp == tempCats(catIdx));
        currentVals = temp_n(temp == tempCats(catIdx));
        temp_p_no = resultsTable.patientNos(fileIdx);
        temp_p_no = temp_p_no{1};
        currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
        ht = text(catIdx,max(currentVals)*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = [label, '_smoker'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(maskFig_n);
    clf;
    temp = resultsTable.usesPatientMask(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    boxplot(temp_n, temp);
    ylabel('Number of particles/m^3/s');
    xlabel('Patient mask');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    tempCats = unique(temp);
    nCats = size(tempCats,1);
    for catIdx = 1:nCats
        currentNevent = nnz(temp == tempCats(catIdx));
        currentVals = temp_n(temp == tempCats(catIdx));
        temp_p_no = resultsTable.patientNos(fileIdx);
        temp_p_no = temp_p_no{1};
        currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
        ht = text(catIdx,max(currentVals)*0.9,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
        set(ht, 'Rotation', 60);
    end
    tempPos = get(gca, 'Position');
    s = 0.9;
    set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
    if saveFigs
        saveFigName = [label, '_mask'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    figure(ageFig_n);
    clf;
    temp = resultsTable.age(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    scatter(temp, temp_n, 'bx', 'LineWidth',1.5);
    ylabel('Number of particles/m^3/s');
    xlabel('Age');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none')
    currentNevent = resultsTable.Nevents(fileIdx);
    currentNpatient = resultsTable.Npatients(fileIdx);
    ht = text(min(temp)*1.05,max(temp_n)*0.9,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
    if saveFigs
        saveFigName = [label, '_age'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    
    figure(bmiFig_n);
    clf;
    temp = resultsTable.bmi(fileIdx);
    temp = temp{1};
    temp_n = resultsTable.n_raw(fileIdx);
    temp_n = temp_n{1};
    scatter(temp, temp_n, 'bx', 'LineWidth',1.5);
    ylabel('Number of particles/m^3/s');
    xlabel('Bmi');
    title(resultsTable.label(fileIdx), 'Interpreter', 'none');
    currentNevent = resultsTable.Nevents(fileIdx);
    currentNpatient = resultsTable.Npatients(fileIdx);
    ht = text(min(temp)*1.05,max(temp_n)*0.9,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
    if saveFigs
        saveFigName = [label, '_bmi'];
        saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
        saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    end
    
    %% Save results table
    save(fullfile(folder,'resultsTable.mat'), 'resultsTable');

end
