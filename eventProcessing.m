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

folder = 'C:\Users\george\OneDrive - The University of Nottingham\SAVE\csv0909_2511';

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
    fileList = {'LowerGI_procedurestarts.csv'};

end

upperGI_nullreffile = 'UpperGI_nullreference.csv';
lowerGI_nullreffile = 'LowerGI_nullreference.csv';

fileList = {upperGI_nullreffile, lowerGI_nullreffile, fileList{1:end}}';

nFiles = size(fileList, 1);
saveFigs = true;
rejectMasks = true;

%% Table for soring vars
varNames = {'label', 'Nevents', 'Npatients', 'mean_n', 'std_n', 'median_n', 'lq_n', 'uq_n',  'mean_v', 'std_v', 'median_v', 'lq_v', 'uq_v', 'mu_mean_n', 'mu_std_n', 'sig_a_n', 'sig_b_n', 'mu_mean_v', 'mu_std_v', 'sig_a_v', 'sig_b_v'};

for k=1:size(varNames,2)
    if k == 1
        varTypes{k} = 'string';
    else
        varTypes{k} = 'double';
    end
end

startFileIdx = 7;
resultsTable = table('Size',[nFiles, size(varNames,2)],'VariableTypes', varTypes, 'VariableNames', varNames);


for fileIdx = 1:nFiles
    file = fileList{fileIdx};
    label = file(1:end-4);
    
    if fileIdx > 2 && fileIdx < startFileIdx
        continue;
    end
    
    if startFileIdx > 1
        if fileIdx == startFileIdx
            load(fullfile(folder,'resultsTable.mat'));
        end
    end
    
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
        roomType_raw = ismember(patientNos,[31,32,33]); %FIX should pull this
        roomType = categorical(roomType_raw+1,1:2,{'Endoscopy room', 'Theatre'});
        
        %UGI route
        ugiRouteTemp_raw = otherVars.RouteOfUGIEndoscopy_oralOrNasal_;
        [ugiRoute, ugiRouteCats] = processVars(otherVars.RouteOfUGIEndoscopy_oralOrNasal_, {'oral', 'nasal', 'N/A'}, {'.*oral.*', '.*mouth.*'}, {'.*nasal.*', '.*nose.*'}, {'.*none.*','.*N/A.*',''});
        
        %Diverticular disease
        diverticularDiseaseTemp_raw = otherVars.ExtensiveDiverticularDisease;
        [diverticularDisease, diverticularDiseaseCats] = processVars(otherVars.ExtensiveDiverticularDisease, {'none', 'mild', 'extensive', 'unknown'}, {'.*none.*', '.*no.*'}, {'.*mild.*', '.*moderate.*'},{'.*extensive.*', '.*yes.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});
        
        %Looping
        loopingTemp_raw = otherVars.DegreeOfLooping;
        [looping, loopingCats] = processVars(otherVars.DegreeOfLooping, {'low', 'medium', 'high', 'unknown/N/A'}, {'.*none.*', '.*low.*', '.*no.*'}, {'.*mild.*', '.*moderate.*', '.*medium.*'},{'.*extensive.*', '.*high.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});
        
        % Discomfort
        discomfortTemp_raw = otherVars.PatientDiscomfort;
        [discomfort, discomfortCats] = processVars(otherVars.PatientDiscomfort, {'low', 'medium', 'high', 'unknown'}, {'.*none.*', '.*low.*', '.*no.*'}, {'.*mild.*', '.*moderate.*', '.*medium.*'},{'.*extensive.*', '.*high.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});
        
        % Hiatus hernia
        hiatusHerniaTemp_raw = otherVars.HiatusHernia;
        if (isnumeric(hiatusHerniaTemp_raw))
            hiatusHernia = categorical(~isnan(hiatusHerniaTemp_raw)+1,1:2,{'no', 'yes'});
            hiatusHerniaCats = categories(hiatusHernia);
        else
            [hiatusHernia, hiatusHerniaCats] = processVars(otherVars.HiatusHernia, {'no', 'yes', 'unknown'}, {'.*none.*', '.*[0].*', '.*no.*'}, {'.*[1-9]?[0-9]*.*', '.*yes.*', '.*massive.*'}, {'.*unknown.*','.*N/A.*',''});
        end
        
        % Suctioning
        suctioningTemp_raw = otherVars.UseOfIntermittentSuctioning;
        [suctioning, suctioningCats] = processVars(otherVars.UseOfIntermittentSuctioning, {'no', 'yes', 'unknown/NA'}, {'.*none.*', '.*no.*'}, {'.*yes.*', '.*y.*'}, {'.*unknown.*','.*N/A.*',''});
       
        a = 1;
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
            realWindowSize = 0;
                   
            for kk=1:size(sample_preInt,2)

                if (kk <= rawdiffwindowSize_after)
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

                    if kk >= rawdiffwindowSize_after
                        break;
                    end
                end

            end

            cumdd = cumdd/realWindowSize;
            sample_int_raw_after_all(:,k) = cumdd;
        end

        %Integrate before
        windowEnd = windowStart;
        rawdiffwindowSize_before = 100;
        windowStart = windowEnd - rawdiffwindowSize_before;

        sampleValid = time > windowStart & time <= (windowStart + rawdiffwindowSize_before + buffer);
        useSmoothedPre = true;
        if useSmoothedPre
            sample_preInt = bg(:,sampleValid,:); %Change to BG to exclude previous evetents
        else
            sample_preInt = raw(:,sampleValid,:); %Change to BG to exclude previous evetents
        end
        sample_int_raw_before_all = zeros(size(raw,1),size(raw,3));

        for k = 1:size(sample_preInt,3)
            cumdt = 0;
            cumdd = zeros(size(sample_preInt,1),1);
            realWindowSize = 0;
            for kk=1:size(sample_preInt,2)

                if (kk <= rawdiffwindowSize_before)
                    cumdt = cumdt + 1;
                end

                if ~all(isnan(sample_preInt(:,kk,k)))
                    if cumdt > avSampleTimes(k)
                        cumdt = avSampleTimes(k); % In case there aren't enough samples before
                    end
                    
                    cumdd = cumdd + sample_preInt(:,kk,k) * cumdt./avSampleTimes(k);

                    realWindowSize = realWindowSize+cumdt;
                    cumdt = 0;

                    if kk >= rawdiffwindowSize_before
                        break;
                    end
                end

            end

            cumdd = cumdd/realWindowSize;
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


        if fileIdx == 1
            curveFitFig = figure;
        else
            figure(curveFitFig);
        end
      
        A_diff = [];
        A_v_diff = [];
        mu_diff = [];
        sigma_diff = [];
        mu_v_diff = [];
        sigma_v_diff = [];
        
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

            [A_est, mu_t, sigma_t, l_t, ~, A_v_est, mu_v_t, sigma_v_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts', 'mu_LB', log(0.1), 'mu_UB', log(20), 'sig_LB', 0.2, 'sig_UB', 50, 'tubeCorrection', correctionVal, 'fullDiameters', diameters_2(:,k));
           
            % Can't use A_est as it artificially inflates larger numbers
            %if A_est == 0
                A_diff(k) = A_t;
                A_v_diff(k) = A_t_v;
            %else
            %    A_diff(k) = A_est;
            %    A_v_diff(k) = A_v_est;
            %end
     
            mu_diff(k) = mu_t;
            sigma_diff(k) = sigma_t;
            mu_v_diff(k) = mu_v_t;
            sigma_v_diff(k) = sigma_v_t;
            
            directFit_v = false;
            if directFit_v
                [A_v_est2, mu_v_t2, sigma_v_t2] = fitAerosolDist(currentDiameters.', (currentData.*currentVols./currentLogBinSizes)','fitType','volume', 'mu_LB', log(0.1), 'mu_UB', log(20), 'sig_LB', 0.2, 'sig_UB', 50, 'tubeCorrection', correctionVal, 'fullDiameters', diameters_2(:,k));
                mu_v_diff(k) = mu_v_t2;
                sigma_v_diff(k) = sigma_v_t2;
            end

        end

        % if clipNegatives
        %     A_fg_after(A_fg_after<0) = 0;
        %     %A_fg_after= A_fg_after(A_fg_after>0);
        % end

        
        if (isUGI_nullref)
            [A_marg_n_mu, A_marg_n_sig] = normfit(A_diff);
            noiseMean_UG = 0;
            noiseStd_UG = A_marg_n_sig;
            
            [A_v_marg_n_mu, A_v_marg_n_sig] = normfit(A_v_diff);
            noiseMean_v_UG = 0;
            noiseStd_v_UG = A_v_marg_n_sig;
            
            noiseMean = noiseMean_UG;
            noiseStd = noiseStd_UG;
            noiseMean_v = noiseMean_v_UG;
            noiseStd_v = noiseStd_v_UG;
            
            reject = false(size(A_diff));
            reject_v = false(size(A_diff));
        elseif isLGI_nullref
            [A_marg_n_mu, A_marg_n_sig] = normfit(A_diff);
            noiseMean_LG = 0;
            noiseStd_LG = A_marg_n_sig;
            
            [A_v_marg_n_mu, A_v_marg_n_sig] = normfit(A_v_diff);
            noiseMean_v_LG = 0;
            noiseStd_v_LG = A_v_marg_n_sig;
            
            
            noiseMean = noiseMean_LG;
            noiseStd = noiseStd_LG;
            noiseMean_v = noiseMean_v_LG;
            noiseStd_v = noiseStd_v_LG;
            
            reject = false(size(A_diff));
            reject_v = false(size(A_diff));
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

                % Fit number
                reject = A_diff < -3*noiseStd; % Reject points that are probably errors
                reject_v = A_v_diff < -3*noiseStd_v; % Reject points that are probably errors
                
                if rejectMasks
                    if strcmpi(label, 'UpperGI_cough') || strcmpi(label, 'UpperGI_procedurestarts') || strcmpi(label, 'UpperGI_procedureends') || strcmpi(label, 'UpperGI_throatspraygiven')
                        reject = reject | (usesPatientMask == 'Mask')';
                        reject_v = reject_v | (usesPatientMask == 'Mask')';
                    end
                    label = [label, '_rejectMasks'];
                end
                
                if nnz(~reject) == 0
                    reject = false(size(reject));
                end

                A_diff_trunc_zeros = A_diff(A_diff>0);
                if isempty(A_diff_trunc_zeros)
                    A_marg_ln_hat = [0,1];
                else
                    A_marg_ln_hat = lognfit(A_diff_trunc_zeros);
                end
                A_diff_trunc = A_diff(~reject);
                objFun = @(x) -1*sum(sumLognormNormpdf(A_diff_trunc, x(1),x(2),0,noiseStd,1e5, exp(log(max(A_diff_trunc))+3*A_marg_ln_hat(2))));
      
                
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
                
                if nnz(~reject_v) == 0
                    reject_v = false(size(reject_v));
                end

                A_v_diff_trunc_zeros = A_v_diff(A_v_diff>0);
                if isempty(A_v_diff_trunc_zeros)
                    A_v_marg_ln_hat = [0,1];
                else
                    A_v_marg_ln_hat = lognfit(A_v_diff_trunc_zeros);
                end
                A_v_diff_trunc = A_v_diff(~reject_v);
                objFun = @(x) -1*sum(sumLognormNormpdf(A_v_diff_trunc, x(1),x(2),0,noiseStd_v, 1e5, exp(log(max(A_v_diff_trunc))+3*A_v_marg_ln_hat(2))));
                
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
                [A_marg_n_mu, A_marg_n_sig] = normfit(A_diff);
            end
        end

        if (~isUGI_nullref && ~isLGI_nullref)
            mu_weights = normcdf(A_diff, noiseMean+3*noiseStd,noiseStd);
            mu_weights = mu_weights/sum(mu_weights)*10;
            [mu_marg_mu, mu_marg_sig] = normfit(mu_diff,[],[],mu_weights);
            sig_marg_hat = gamfit(sigma_diff,[],[],mu_weights);
            
            mu_v_weights = normcdf(A_v_diff, noiseMean_v+3*noiseStd_v,noiseStd_v);
            mu_v_weights = mu_v_weights/sum(mu_v_weights)*10;
            [mu_v_marg_mu, mu_v_marg_sig] = normfit(mu_v_diff,[],[],mu_v_weights*100);
            sig_v_marg_hat = gamfit(sigma_v_diff,[],[],mu_v_weights*100);
        else
            [mu_marg_mu, mu_marg_sig] = normfit(mu_diff);
            sig_marg_hat = gamfit(sigma_diff);
            
            [mu_v_marg_mu, mu_v_marg_sig] = normfit(mu_v_diff);
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
        scatter(exp(mu_v_diff),zeros(size(mu_v_diff)),'rx','LineWidth',2);
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
        scatter(A_v_diff,zeros(size(A_v_diff)),'kx','LineWidth',2);
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
        scatter(A_v_diff, exp(mu_v_diff),'gx','LineWidth',2);
        xlim([min(A_v_plot), max(A_v_plot)]);
        ylim([min(exp(mu_plot)), max(exp(mu_plot))]);
        xlabel('vol particles/m^3/s');
        ylabel('particle diameter (\mu m)');
        title('Joint distribution \mu vs #');

        subplot(3,2,4);
        scatter(exp(mu_v_diff),sigma_v_diff,'gx','LineWidth',2);
        xlim([min(exp(mu_plot)), max(exp(mu_plot))]);
        ylim([min(sig_plot), max(sig_plot)]);
        xlabel('particle diameter (\mu m)');
        ylabel('sigma');
        title('Joint distribution \sigma vs \mu');

        subplot(3,2,6);
        scatter(A_v_diff,sigma_v_diff,'gx','LineWidth',2);
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
    nEvents = size(indices(~reject),1);
    nPatients = size(unique(patientNos(~reject)),1);
    
    nEvents_v = size(indices(~reject_v),1);
    nPatients_v = size(unique(patientNos(~reject_v)),1);
    
    
    resultsTable.label(fileIdx) = label;
    resultsTable.Nevents(fileIdx) = nEvents;
    resultsTable.Npatients(fileIdx) = nPatients;
    resultsTable.Nevents_v(fileIdx) = nEvents_v;
    resultsTable.Npatients_v(fileIdx) = nPatients_v;
    
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
        else %Delete
            resultsTable.mean_n(fileIdx) = A_marg_n_mu;
            resultsTable.std_n(fileIdx) = A_marg_n_sigma;
            resultsTable.mean_v(fileIdx) = A_v_marg_n_mu;
            resultsTable.std_v(fileIdx) = A_v_marg_n_sigma;
        end
    end
    
    resultsTable.median_n(fileIdx) = median(A_diff(~reject));
    resultsTable.lq_n(fileIdx) = quantile(A_diff(~reject), 0.25);
    resultsTable.uq_n(fileIdx) = quantile(A_diff(~reject), 0.75);
      
    resultsTable.median_v(fileIdx) = median(A_v_diff(~reject_v));
    resultsTable.lq_v(fileIdx) = quantile(A_v_diff(~reject_v), 0.25);
    resultsTable.uq_v(fileIdx) = quantile(A_v_diff(~reject_v), 0.75);
    
    resultsTable.mu_mean_n(fileIdx) = mu_marg_mu;
    resultsTable.mu_std_n(fileIdx) = mu_marg_sig;  
    resultsTable.mu_mean_v(fileIdx) = mu_v_marg_mu;
    resultsTable.mu_std_v(fileIdx) = mu_v_marg_sig;
      
    resultsTable.sig_a_n(fileIdx) = sig_marg_hat(1);
    resultsTable.sig_b_n(fileIdx) = sig_marg_hat(2);
    resultsTable.sig_a_v(fileIdx) = sig_v_marg_hat(1);
    resultsTable.sig_b_v(fileIdx) = sig_v_marg_hat(2);
    
    resultsTable.n_raw(fileIdx) = {A_diff(~reject)}; %Truncate
    resultsTable.v_raw(fileIdx) = {A_v_diff(~reject_v)};
    resultsTable.mu_raw(fileIdx) = {mu_diff(~reject)};
    resultsTable.mu_v_raw(fileIdx) = {mu_v_diff(~reject_v)};
    
    resultsTable.age(fileIdx) = {otherVars.Age(~reject)};
    resultsTable.bmi(fileIdx) = {otherVars.BMI(~reject)};
    resultsTable.sedation(fileIdx) = {sedation(~reject)};
    resultsTable.analTone(fileIdx) = {analTone(~reject)};
    resultsTable.useOfCO2OrWater(fileIdx) = {useOfCO2orWater(~reject)};
    resultsTable.isSmoker(fileIdx) = {isSmoker(~reject)};
    resultsTable.usesPatientMask(fileIdx) = {usesPatientMask(~reject)};
    resultsTable.patientNos(fileIdx) = {patientNos(~reject)};
    resultsTable.roomType(fileIdx) = {roomType(~reject)};
    resultsTable.ugiRoute(fileIdx) = {ugiRoute(~reject)};
    resultsTable.diverticularDisease(fileIdx) = {diverticularDisease(~reject)};
    resultsTable.looping(fileIdx) = {looping(~reject)};
    resultsTable.discomfort(fileIdx) = {discomfort(~reject)};
    resultsTable.hiatusHernia(fileIdx) = {hiatusHernia(~reject)};
    resultsTable.suctioning(fileIdx) = {suctioning(~reject)};
    
    resultsTable.age_v(fileIdx) = {otherVars.Age(~reject_v)};
    resultsTable.bmi_v(fileIdx) = {otherVars.BMI(~reject_v)};
    resultsTable.sedation_v(fileIdx) = {sedation(~reject_v)};
    resultsTable.analTone_v(fileIdx) = {analTone(~reject_v)};
    resultsTable.useOfCO2OrWater_v(fileIdx) = {useOfCO2orWater(~reject_v)};
    resultsTable.isSmoker_v(fileIdx) = {isSmoker(~reject_v)};
    resultsTable.usesPatientMask_v(fileIdx) = {usesPatientMask(~reject_v)};
    resultsTable.patientNos_v(fileIdx) = {patientNos(~reject_v)};
    resultsTable.roomType_v(fileIdx) = {roomType(~reject_v)};
    resultsTable.ugiRoute_v(fileIdx) = {ugiRoute(~reject_v)};
    resultsTable.diverticularDisease_v(fileIdx) = {diverticularDisease(~reject_v)};
    resultsTable.looping_v(fileIdx) = {looping(~reject_v)};
    resultsTable.discomfort_v(fileIdx) = {discomfort(~reject_v)};
    resultsTable.hiatusHernia_v(fileIdx) = {hiatusHernia(~reject_v)};
    resultsTable.suctioning_v(fileIdx) = {suctioning(~reject_v)};
    
    resultsTable.samples{fileIdx} = []; %Delete
    
    % Compute pvals
    for tempIdx = 1:fileIdx-1
        event1 = resultsTable.n_raw{tempIdx};
        event1 = event1(:);
        event2 = resultsTable.n_raw{fileIdx};
        event2 = event2(:);
        
        %FIX what to do about first distribution
        if tempIdx == 1 || tempIdx == 2
            mean1 = -1.99;
            std1 = 0.1;
        else
            mean1 = resultsTable.mean_n(tempIdx);
            std1 = resultsTable.std_n(tempIdx);
        end
        
        if fileIdx == 1 || fileIdx == 2
            mean2 = -1.99;
            std2 = 0.1;
        else
            mean2 = resultsTable.mean_n(fileIdx);
            std2 = resultsTable.std_n(fileIdx);
        end
        
        computeEventPvals = true;    
        if computeEventPvals
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [pMu, pSig, samples1out, samples2out] = computeSignificance(event1, event2, noiseMean, noiseStd, mean1, std1, mean2, std2);
        else
            pMu = 0.5;
            pSig = 0.5;
        end
            
        pMu = min([pMu, 1-pMu]);
        pSig = min([pSig, 1-pSig]);
     
        pMuTable(fileIdx, tempIdx) = pMu;
        pSigTable(fileIdx, tempIdx) = pSig;
        
        
        % Now volume
        event1 = resultsTable.v_raw{tempIdx};
        event1 = event1(:);
        event2 = resultsTable.v_raw{fileIdx};
        event2 = event2(:);
        
        
        if tempIdx == 1 || tempIdx == 2
            mean1 = -32.9;
            std1 = 0.1;
        else
            mean1 = resultsTable.mean_v(tempIdx);
            std1 = resultsTable.std_v(tempIdx);
        end
        
        if fileIdx == 1 || fileIdx == 2
            mean2 = -32.9;
            std2 = 0.1;
        else
            mean2 = resultsTable.mean_v(fileIdx);
            std2 = resultsTable.std_v(fileIdx);
        end
        
        computeEventPvals_v = true;    
        if computeEventPvals_v
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [pMu_v, pSig_v] = computeSignificance(event1, event2, noiseMean_v, noiseStd_v, mean1, std1, mean2, std2, 'muMinIn', -34);
        else
            pMu_v = 0.5;
            pSig_v = 0.5;
        end
            
        pMu_v = min([pMu_v, 1-pMu_v]);
        pSig_v = min([pSig_v, 1-pSig_v]);
     
        pMuTable_v(fileIdx, tempIdx) = pMu_v;
        pSigTable_v(fileIdx, tempIdx) = pSig_v;
        
    end
    
    data_box = [];
    data_v_box = [];
    mu_box = [];
    mu_v_box = [];
    cats_box = [];
    cats_v_box = [];
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
        
        current_cats_v = repmat(resultsTable.label(dataSetIdx),size(currentRaw_v,1),1);
        current_cats_v = categorical(current_cats_v);
        
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
        cats_v_box = [cats_v_box; current_cats_v;];
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
        
        
        dataCompFigs(1).fig = boxFig_n;
        dataCompFigs(1).plotType = 'boxplot';
        dataCompFigs(1).varNames = 'data_box';
        dataCompFigs(1).catNames = 'cats_box';
        dataCompFigs(1).plotPvals = true;
        dataCompFigs(1).pValVarMu = 'pMuTable';
        dataCompFigs(1).pValVarSig = 'pSigTable';
        dataCompFigs(1).plotNonLogYticks = false;
        dataCompFigs(1).rawData = 'n_raw';
        dataCompFigs(1).meanVar = 'mean_n';
        dataCompFigs(1).ylabel = 'Number of particles/m^3/s';
        dataCompFigs(1).saveName = ['particle_number_', dataCompFigs(1).plotType];
        
        dataCompFigs(2).fig = violinFig_n;
        dataCompFigs(2).plotType = 'violinplot';
        dataCompFigs(2).varNames = 'data_box';
        dataCompFigs(2).catNames = 'cats_box';
        dataCompFigs(2).plotPvals = true;
        dataCompFigs(2).pValVarMu = 'pMuTable';
        dataCompFigs(2).pValVarSig = 'pSigTable';
        dataCompFigs(2).plotNonLogYticks = false;
        dataCompFigs(2).rawData = 'n_raw';
        dataCompFigs(2).meanVar = 'mean_n';
        dataCompFigs(2).ylabel = 'Number of particles/m^3/s';
        dataCompFigs(2).saveName = ['particle_number_', dataCompFigs(2).plotType];
        
        dataCompFigs(3).fig = boxFig_v;
        dataCompFigs(3).plotType = 'boxplot';
        dataCompFigs(3).varNames = 'data_v_box';
        dataCompFigs(3).catNames = 'cats_v_box';
        dataCompFigs(3).plotPvals = true;
        dataCompFigs(3).pValVarMu = 'pMuTable_v';
        dataCompFigs(3).pValVarSig = 'pSigTable_v';
        dataCompFigs(3).plotNonLogYticks = false;
        dataCompFigs(3).rawData = 'v_raw';
        dataCompFigs(3).meanVar = 'mean_v';
        dataCompFigs(3).ylabel = 'Volume of particles (m^3) /m^3/s';
        dataCompFigs(3).saveName = ['particle_volume_', dataCompFigs(3).plotType];
        
        dataCompFigs(4).fig = violinFig_v;
        dataCompFigs(4).plotType = 'violinplot';
        dataCompFigs(4).varNames = 'data_v_box';
        dataCompFigs(4).catNames = 'cats_v_box';
        dataCompFigs(4).plotPvals = true;
        dataCompFigs(4).pValVarMu = 'pMuTable_v';
        dataCompFigs(4).pValVarSig = 'pSigTable_v';
        dataCompFigs(4).plotNonLogYticks = false;
        dataCompFigs(4).rawData = 'v_raw';
        dataCompFigs(4).meanVar = 'mean_v';
        dataCompFigs(4).ylabel = 'Volume of particles (m^3) /m^3/s';
        dataCompFigs(4).saveName = ['particle_volume_', dataCompFigs(4).plotType];
        
        dataCompFigs(5).fig = boxFig_mu_n;
        dataCompFigs(5).plotType = 'boxplot';
        dataCompFigs(5).varNames = 'mu_box';
        dataCompFigs(5).catNames = 'cats_box';
        dataCompFigs(5).plotPvals = false;
        dataCompFigs(5).pValVarMu = '';
        dataCompFigs(5).pValVarSig = '';
        dataCompFigs(5).plotNonLogYticks = true;
        dataCompFigs(5).rawData = 'mu_raw';
        dataCompFigs(5).meanVar = '';
        dataCompFigs(5).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(5).saveName = ['diameter_number_', dataCompFigs(5).plotType];
        
        dataCompFigs(6).fig = violinFig_mu_n;
        dataCompFigs(6).plotType = 'violinplot';
        dataCompFigs(6).varNames = 'mu_box';
        dataCompFigs(6).catNames = 'cats_box';
        dataCompFigs(6).plotPvals = false;
        dataCompFigs(6).pValVarMu = '';
        dataCompFigs(6).pValVarSig = '';
        dataCompFigs(6).plotNonLogYticks = true;
        dataCompFigs(6).rawData = 'mu_raw';
        dataCompFigs(6).meanVar = '';
        dataCompFigs(6).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(6).saveName = ['diameter_number_', dataCompFigs(6).plotType];
        
        dataCompFigs(7).fig = boxFig_mu_v;
        dataCompFigs(7).plotType = 'boxplot';
        dataCompFigs(7).varNames = 'mu_v_box';
        dataCompFigs(7).catNames = 'cats_v_box';
        dataCompFigs(7).plotPvals = false;
        dataCompFigs(7).pValVarMu = '';
        dataCompFigs(7).pValVarSig = '';
        dataCompFigs(7).plotNonLogYticks = true;
        dataCompFigs(7).rawData = 'mu_v_raw';
        dataCompFigs(7).meanVar = '';
        dataCompFigs(7).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(7).saveName = ['diameter_volume_', dataCompFigs(7).plotType];
        
        dataCompFigs(8).fig = violinFig_mu_v;
        dataCompFigs(8).plotType = 'violinplot';
        dataCompFigs(8).varNames = 'mu_v_box';
        dataCompFigs(8).catNames = 'cats_v_box';
        dataCompFigs(8).plotPvals = false;
        dataCompFigs(8).pValVarMu = '';
        dataCompFigs(8).pValVarSig = '';
        dataCompFigs(8).plotNonLogYticks = true;
        dataCompFigs(8).rawData = 'mu_v_raw';
        dataCompFigs(8).meanVar = '';
        dataCompFigs(8).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(8).saveName = ['diameter_volume_', dataCompFigs(8).plotType];
    end
    
    for k=1:size(dataCompFigs,2)
        figure(dataCompFigs(k).fig);
        clf;
        eval([dataCompFigs(k).plotType, '(', dataCompFigs(k).varNames, ',', dataCompFigs(k).catNames, ');']);
        tempCats_ = eval(['unique(', dataCompFigs(k).catNames, ')']);
        set(gca,'TickLabelInterpreter', 'none');
        xtickangle(60);
        hold on;
        
        if (k <= 4)
            meanVals = exp(eval(['resultsTable.', dataCompFigs(k).meanVar, '(1:fileIdx)']));
            scatter(1:fileIdx, meanVals, 150,'kx', 'LineWidth', 3); 
        end
        
       	if dataCompFigs(k).plotNonLogYticks
            y_ticks = yticklabels;
            for yidx = 1:size(y_ticks,1)
                temp = str2num(y_ticks{yidx});
                temp = exp(temp);
                temp = sprintf('%2.2f',temp);
                y_ticks{yidx} = temp;
            end
            yticklabels(y_ticks);
        end
        ylabel(dataCompFigs(k).ylabel);
        for tempIdx = 1:fileIdx
            currentNevent = resultsTable.Nevents(tempIdx);
            currentNpatient = resultsTable.Npatients(tempIdx);
            ht = text(tempIdx,max(data_box(cats_box == tempCats_(tempIdx)))*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
            set(ht,'Rotation',60);

            if dataCompFigs(k).plotPvals
                if tempIdx < fileIdx
                    for fileIdx2 = tempIdx+1:fileIdx
                        tempD1 = eval(['resultsTable.', dataCompFigs(k).rawData, '{tempIdx};']);
                        tempD1 = tempD1(:);
                        tempD2 = eval(['resultsTable.', dataCompFigs(k).rawData, '{fileIdx2};']);
                        tempD2 = tempD2(:);

                        if (fileIdx > 1) && (tempIdx < fileIdx)
                            pValMu = eval([dataCompFigs(k).pValVarMu, '(fileIdx2, tempIdx);']);
                            pValSig = eval([dataCompFigs(k).pValVarSig, '(fileIdx2, tempIdx);']);
                            computeAndPlotPvals(tempD1,tempD2,[],[],tempIdx, fileIdx2,'pMeanIn', pValMu, 'pSigIn', pValSig);
                        end
                    end
                end
            end
        end   
        tempPos = get(gca, 'Position');
        s = 0.9;
        set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
        if saveFigs
            saveFigName = dataCompFigs(k).saveName;
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
    end
    
    %% Now plot variables
    if (fileIdx == 1)
        sedationFig_n = figure;
        analToneFig_n = figure;
        co2vWaterFig_n = figure;
        smokerFig_n = figure;
        maskFig_n = figure;
        roomTypeFig_n = figure;
        ugiRouteFig_n = figure;
        divertDiseaseFig_n = figure;
        loopingFig_n = figure;
        discomfortFig_n = figure;
        hiatusFig_n = figure;
        suctioningFig_n = figure;
        ageFig_n = figure;
        bmiFig_n = figure;
        
        sedationFig_v = figure;
        analToneFig_v = figure;
        co2vWaterFig_v = figure;
        smokerFig_v = figure;
        maskFig_v = figure;
        roomTypeFig_v = figure;
        ugiRouteFig_v = figure;
        divertDiseaseFig_v = figure;
        loopingFig_v = figure;
        discomfortFig_v = figure;
        hiatusFig_v = figure;
        suctioningFig_v = figure;
        ageFig_v = figure;
        bmiFig_v = figure;
        
        varFigs(1).fig = sedationFig_n;
        varFigs(1).varname = 'sedation';
        varFigs(1).rawvar = 'n_raw';
        varFigs(1).type = 'discrete';
        varFigs(2).fig = analToneFig_n;
        varFigs(2).varname = 'analTone';
        varFigs(2).rawvar = 'n_raw';
        varFigs(2).type = 'discrete';
        varFigs(3).fig = co2vWaterFig_n;
        varFigs(3).varname = 'useOfCO2OrWater';
        varFigs(3).rawvar = 'n_raw';
        varFigs(3).type = 'discrete';
        varFigs(4).fig = smokerFig_n;
        varFigs(4).varname = 'isSmoker';
        varFigs(4).rawvar = 'n_raw';
        varFigs(4).type = 'discrete';
        varFigs(5).fig = maskFig_n;
        varFigs(5).varname = 'usesPatientMask';
        varFigs(5).rawvar = 'n_raw';
        varFigs(5).type = 'discrete';
        varFigs(6).fig = roomTypeFig_n;
        varFigs(6).varname = 'roomType';
        varFigs(6).rawvar = 'n_raw';
        varFigs(6).type = 'discrete';
        varFigs(7).fig = ugiRouteFig_n;
        varFigs(7).varname = 'ugiRoute';
        varFigs(7).type = 'discrete';
        varFigs(7).rawvar = 'n_raw';
        varFigs(8).fig = divertDiseaseFig_n;
        varFigs(8).varname = 'diverticularDisease';
        varFigs(8).rawvar = 'n_raw';
        varFigs(8).type = 'discrete';
        varFigs(9).fig = loopingFig_n;
        varFigs(9).varname = 'looping';
        varFigs(9).rawvar = 'n_raw';
        varFigs(9).type = 'discrete';
        varFigs(10).fig = discomfortFig_n;
        varFigs(10).varname = 'discomfort';
        varFigs(10).rawvar = 'n_raw';
        varFigs(10).type = 'discrete';
        varFigs(11).fig = hiatusFig_n;
        varFigs(11).varname = 'hiatusHernia';
        varFigs(11).rawvar = 'n_raw';
        varFigs(11).type = 'discrete';
        varFigs(12).fig = suctioningFig_n;
        varFigs(12).varname = 'suctioning';
        varFigs(12).rawvar = 'n_raw';
        varFigs(12).type = 'discrete';
        varFigs(13).fig = ageFig_n;
        varFigs(13).varname = 'age';
        varFigs(13).rawvar = 'n_raw';
        varFigs(13).type = 'continuous';
        varFigs(14).fig = bmiFig_n;
        varFigs(14).varname = 'bmi';
        varFigs(14).rawvar = 'n_raw';
        varFigs(14).type = 'continuous';
        
        %repeat for volume
        varFigs(15).fig = sedationFig_v;
        varFigs(15).varname = 'sedation_v';
        varFigs(15).rawvar = 'v_raw';
        varFigs(15).type = 'discrete';
        varFigs(16).fig = analToneFig_v;
        varFigs(16).varname = 'analTone_v';
        varFigs(16).rawvar = 'v_raw';
        varFigs(16).type = 'discrete';
        varFigs(17).fig = co2vWaterFig_v;
        varFigs(17).varname = 'useOfCO2OrWater_v';
        varFigs(17).rawvar = 'v_raw';
        varFigs(17).type = 'discrete';
        varFigs(18).fig = smokerFig_v;
        varFigs(18).varname = 'isSmoker_v';
        varFigs(18).rawvar = 'v_raw';
        varFigs(18).type = 'discrete';
        varFigs(19).fig = maskFig_v;
        varFigs(19).varname = 'usesPatientMask_v';
        varFigs(19).rawvar = 'v_raw';
        varFigs(19).type = 'discrete';
        varFigs(20).fig = roomTypeFig_v;
        varFigs(20).varname = 'roomType_v';
        varFigs(20).rawvar = 'v_raw';
        varFigs(20).type = 'discrete';
        varFigs(21).fig = ugiRouteFig_v;
        varFigs(21).varname = 'ugiRoute_v';
        varFigs(21).rawvar = 'v_raw';
        varFigs(21).type = 'discrete';
        varFigs(22).fig = divertDiseaseFig_v;
        varFigs(22).varname = 'diverticularDisease_v';
        varFigs(22).rawvar = 'v_raw';
        varFigs(22).type = 'discrete';
        varFigs(23).fig = loopingFig_v;
        varFigs(23).varname = 'looping_v';
        varFigs(23).rawvar = 'v_raw';
        varFigs(23).type = 'discrete';
        varFigs(24).fig = discomfortFig_v;
        varFigs(24).varname = 'discomfort_v';
        varFigs(24).rawvar = 'v_raw';
        varFigs(24).type = 'discrete';
        varFigs(25).fig = hiatusFig_v;
        varFigs(25).varname = 'hiatusHernia_v';
        varFigs(25).rawvar = 'v_raw';
        varFigs(25).type = 'discrete';
        varFigs(26).fig = suctioningFig_v;
        varFigs(26).varname = 'suctioning_v';
        varFigs(26).rawvar = 'v_raw';
        varFigs(26).type = 'discrete';
        varFigs(27).fig = ageFig_v;
        varFigs(27).varname = 'age_v';
        varFigs(27).rawvar = 'v_raw';
        varFigs(27).type = 'continuous';
        varFigs(28).fig = bmiFig_v;
        varFigs(28).varname = 'bmi_v';
        varFigs(28).rawvar = 'v_raw';
        varFigs(28).type = 'continuous';
    end
    
    for k=1:size(varFigs,2)
        figure(varFigs(k).fig);
        clf;
        eval(['temp = resultsTable.', varFigs(k).varname, '(fileIdx);']);
        temp = temp{1};
        temp_n = eval(['resultsTable.', varFigs(k).rawvar, '(fileIdx)']);
        temp_n = temp_n{1};
        if strcmpi(varFigs(k).type, 'discrete')
            %boxplot(temp_n, temp);
            violinplot(temp_n, temp);
            tempCats = unique(temp);
            nCats = size(tempCats,1);
            for catIdx = 1:nCats
                currentNevent = nnz(temp == tempCats(catIdx));
                currentVals = temp_n(temp == tempCats(catIdx));
                if strcmpi(varFigs(k).rawvar, 'n_raw')
                	temp_p_no = resultsTable.patientNos(fileIdx);
                elseif strcmpi(varFigs(k).rawvar, 'v_raw')
                	temp_p_no = resultsTable.patientNos_v(fileIdx);
                end
                temp_p_no = temp_p_no{1};
                currentNpatient = nnz(unique(temp_p_no(temp == tempCats(catIdx))));
                ht = text(catIdx,max(currentVals)*1.05,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);
                set(ht, 'Rotation', 60);
            end
            if nCats > 1
                for cat1Idx = 1:nCats-1
                     for cat2Idx = cat1Idx+1:nCats                 
                         val1 = temp == tempCats(cat1Idx);
                         val2 = temp == tempCats(cat2Idx);

                         d1 = temp_n(val1);
                         d2 = temp_n(val2);

                         if fileIdx > 2 % && strcmpi(varFigs(k).varname, 'usesPatientMask')
                            if strcmpi(varFigs(k).rawvar, 'n_raw')
                                computeAndPlotPvals(d1,d2,noiseMean,noiseStd,cat1Idx,cat2Idx, 'pThresh', 0.5);
                            elseif strcmpi(varFigs(k).rawvar, 'v_raw')
                                computeAndPlotPvals(d1,d2,noiseMean,noiseStd,cat1Idx,cat2Idx, 'pThresh', 0.5, 'muMinIn', -34);
                            end
                         end
                     end
                end
            end
            tempPos = get(gca, 'Position');
            s = 0.9;
            set(gca,'Position', [tempPos(1)+tempPos(3)*(1-s)/2, tempPos(2)+tempPos(4)*(1-s)/2,tempPos(3)*s, tempPos(4)*s]);
        else
            scatter(temp, temp_n, 'bx', 'LineWidth',1.5);
            title(resultsTable.label(fileIdx), 'Interpreter', 'none')
            if strcmpi(varFigs(k).rawvar, 'n_raw')
                currentNevent = resultsTable.Nevents(fileIdx);
                currentNpatient = resultsTable.Npatients(fileIdx);
            elseif strcmpi(varFigs(k).rawvar, 'v_raw')
                currentNevent = resultsTable.Nevents_v(fileIdx);
                currentNpatient = resultsTable.Npatients_v(fileIdx);
            end
            valid = ~isnan(temp);
            temp = temp(valid);
            temp_n = temp_n(valid);
            correctionVal = nnz(~valid);
            currentNpatient = currentNpatient - correctionVal;
            currentNevent = currentNevent - correctionVal;
            ht = text(min(temp)*1.05,max(temp_n)*0.9,['N_{pat} = ', num2str(currentNpatient), ', N_{ev} = ', num2str(currentNevent)]);

            if (size(temp,1) > 1)
                newx = linspace(min(temp),max(temp),100);
                [fitresult, gof, ~] = fit(temp(:),temp_n(:),'poly1');
                yfit = feval(fitresult,newx);

                if (size(temp,1) > 2)
                    p21 = predint(fitresult,newx,0.95,'functional','off');
                end
                hold on;
                plot(newx, yfit, 'k');

                if (size(temp,1) > 2)
                    plot(newx, p21, 'm--');
                end
                text(min(newx)*1.1, max(temp_n)*0.8,['r = ', num2str(gof.rsquare)]);
                hold off;
            end
        end
        
        if strcmpi(varFigs(k).rawvar, 'n_raw')
            ylabel('Number of particles/m^3/s');
        elseif strcmpi(varFigs(k).rawvar, 'v_raw')
            ylabel('Volume of particles m^3/m^3/s');
        end
        xlabel(varFigs(k).varname);
        title(resultsTable.label(fileIdx), 'Interpreter', 'none');
        
        if saveFigs
            saveFigName = [label, '_', varFigs(k).varname];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
    end
    
    %%
    
    if (fileIdx == 1)
        dTree_fig = figure;
        forest_fig = figure;
    end
        
    if (fileIdx > 2)
    
        % Make a table to fit
        n_raw_tab = resultsTable.n_raw{fileIdx}';
        v_raw_tab = resultsTable.v_raw{fileIdx}';
        mu_raw_tab = resultsTable.mu_raw{fileIdx}';
        mu_v_raw_tab = resultsTable.mu_v_raw{fileIdx}';
        age_tab = resultsTable.age{fileIdx};
        bmi_tab = resultsTable.bmi{fileIdx};
        sedation_tab = resultsTable.sedation{fileIdx};
        analTone_tab = resultsTable.analTone{fileIdx};
        co2water_tab = resultsTable.useOfCO2OrWater{fileIdx};
        isSmoker_tab = resultsTable.isSmoker{fileIdx};
        patientMask_tab = resultsTable.usesPatientMask{fileIdx};
        %patientNos_tab = resultsTable.patientNos{fileIdx};
        roomType_tab = resultsTable.roomType{fileIdx};
        ugiRoute_tab = resultsTable.ugiRoute{fileIdx};
        divertDisease_tab = resultsTable.diverticularDisease{fileIdx};
        looping_tab = resultsTable.looping{fileIdx};
        discomfort_tab = resultsTable.discomfort{fileIdx};
        hiatusHernia_tab = resultsTable.hiatusHernia{fileIdx};
        suctioning_tab = resultsTable.suctioning{fileIdx};

        tempTab = table(age_tab, bmi_tab, sedation_tab, analTone_tab, co2water_tab, isSmoker_tab, patientMask_tab, roomType_tab, ugiRoute_tab, divertDisease_tab, looping_tab, discomfort_tab, hiatusHernia_tab, suctioning_tab);

        tempTree = fitrtree(tempTab, n_raw_tab, 'CategoricalPredictors',[3:size(table,2)],'Surrogate','on');
        before = findall(groot,'Type','figure'); % Find all figures
        view(tempTree,'Mode','graph');
        after = findall(groot,'Type','figure');
        h = setdiff(after,before); % Get the figure handle of the tree viewer
        if saveFigs
            saveFigName = [label, '_regression_tree'];
            saveas(h,fullfile(folder,[saveFigName, '.fig']));
            saveas(h,fullfile(folder,[saveFigName, '.png']));
        end
        close(h);
        
        tempForest = TreeBagger(1000, tempTab, n_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on');

        %view(tempForest.Trees{1},'Mode','graph')
        imp = tempForest.OOBPermutedPredictorDeltaError;
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,2,1);
        bar(imp);
        title('Variable importance');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,2,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        if saveFigs
            saveFigName = [label, '_random_forest'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
        
        %% Now for volumes
        age_tab_v = resultsTable.age_v{fileIdx};
        bmi_tab_v = resultsTable.bmi_v{fileIdx};
        sedation_tab_v = resultsTable.sedation_v{fileIdx};
        analTone_tab_v = resultsTable.analTone_v{fileIdx};
        co2water_tab_v = resultsTable.useOfCO2OrWater_v{fileIdx};
        isSmoker_tab_v = resultsTable.isSmoker_v{fileIdx};
        patientMask_tab_v = resultsTable.usesPatientMask_v{fileIdx};
        %patientNos_tab_v = resultsTable.patientNos_v{fileIdx};
        roomType_tab_v = resultsTable.roomType_v{fileIdx};
        ugiRoute_tab_v = resultsTable.ugiRoute_v{fileIdx};
        divertDisease_tab_v = resultsTable.diverticularDisease_v{fileIdx};
        looping_tab_v = resultsTable.looping_v{fileIdx};
        discomfort_tab_v = resultsTable.discomfort_v{fileIdx};
        hiatusHernia_tab_v = resultsTable.hiatusHernia_v{fileIdx};
        suctioning_tab_v = resultsTable.suctioning_v{fileIdx};

        tempTab_v = table(age_tab_v, bmi_tab_v, sedation_tab_v, analTone_tab_v, co2water_tab_v, isSmoker_tab_v, patientMask_tab_v, roomType_tab_v, ugiRoute_tab_v, divertDisease_tab_v, looping_tab_v, discomfort_tab_v, hiatusHernia_tab_v, suctioning_tab_v);

        tempTree = fitrtree(tempTab_v, v_raw_tab, 'CategoricalPredictors',[3:size(table,2)],'Surrogate','on');
        before = findall(groot,'Type','figure'); % Find all figures
        view(tempTree,'Mode','graph');
        after = findall(groot,'Type','figure');
        h = setdiff(after,before); % Get the figure handle of the tree viewer
        if saveFigs
            saveFigName = [label, '_regression_tree_v'];
            saveas(h,fullfile(folder,[saveFigName, '.fig']));
            saveas(h,fullfile(folder,[saveFigName, '.png']));
        end
        close(h);
        
        tempForest = TreeBagger(1000, tempTab_v, v_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on');

        %view(tempForest.Trees{1},'Mode','graph')
        imp = tempForest.OOBPermutedPredictorDeltaError;
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,2,1);
        bar(imp);
        title('Variable importance vol.');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,2,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        if saveFigs
            saveFigName = [label, '_random_forest_v'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
    end

    %% Save results table
    if fileIdx >= startFileIdx
        if exist('pMuTable')
            save(fullfile(folder,'resultsTable.mat'), 'resultsTable', 'pMuTable', 'pSigTable', 'pMuTable_v', 'pSigTable_v');
        else
            save(fullfile(folder,'resultsTable.mat'), 'resultsTable');
        end
    end
    
    %%
    useFGBGapproach = false;
    if useFGBGapproach%% fit FG after the event
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
    end

    %% STAN
    useSTAN = false;
    if useSTAN
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
    end

end
