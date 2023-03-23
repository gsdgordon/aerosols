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

useCSVs = false;
saveFigs = true;
rejectMasks = true;
rejectERCP = true;
rejectLF = true;
rejectTheatre = true;
LFonly = false;
rejectNasalAbandoned = true;
cytospongeOnly = true;
rejectCytosponge = false;

limitSize = true;
lt = false;
sizeLim = 5;
analyseVariables = true;

computeEventPvals = true; 
computeEventPvals_v = false;    
computeEventPvals_mu = false;    
computeEventPvals_sig = false;
computeVarPvals = false;

username = getenv('username');

if useCSVs
    folder = ['C:\Users\', username, '\OneDrive - The University of Nottingham\SAVE\csv0909_2511'];

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
else
    % This loads all events and times
    folder = sprintf('%s_results', datestr(now,'yyyy-mm-dd'));
    %folder = '01-07-2021_results';
    
    dataDir = ['C:\Users\', username, '\OneDrive - The University of Nottingham\SAVE\'];
    
    folder = [dataDir, folder];
    if ~exist(folder)
        mkdir(folder);
    end
    
%     eventList = {'Upper GI_Null reference',...
%                  'Lower GI_Null reference',...
%                 'Upper GI_Deep breaths',...
%                 'Upper GI_Forced cough',...
%                 'Upper GI_Speaking',...
%                 'Upper GI_Procedure starts',...
%                 'Upper GI_Procedure ends',...
%                 'Upper GI_Procedure starts nasal',...
%                 'Upper GI_Procedure ends nasal',...
%                 'Upper GI_Cough',...
%                 'Upper GI_Cough during nasal',...
%                 'Upper GI_Fundal retroflexion',...
%                 'Upper GI_Nasal spray given',...
%                 'Upper GI_Throat spray given',...
%                 'Lower GI_Procedure starts',...
%                 'Lower GI_Procedure ends',...
%                 'Lower GI_Abdominal pressure',...
%                 'Lower GI_Biopsy sampling',...
%                 'Lower GI_Caecal intubation',...
%                 'Lower GI_Catheter in',...
%                 'Lower GI_Catheter out',...
%                 'Lower GI_Position changes',...
%                 'Lower GI_Diathermy pad used',...
%                 'Lower GI_Rectal insufflation/retroflexion'};
% % 
%     eventList = {'Upper GI_Null reference',...
%                  'Lower GI_Null reference',...
%                 'Upper GI_Forced cough',...
%                 'Upper GI_Deep breaths',...
%                 'Upper GI_Speaking',...
%                 'Upper GI_Procedure starts',...
%                 'Upper GI_Procedure ends',...
%                 'Upper GI_Procedure starts nasal',...
%                 'Upper GI_Procedure ends nasal',...
%                 'Upper GI_Cough',...
%                 'Upper GI_Cough during nasal',...
%                 'Upper GI_Fundal retroflexion',...
%                 'Upper GI_Nasal spray given',...
%                 'Upper GI_Throat spray given',...
%                 'Upper GI_Biopsy sampling'};
            
%      eventList = {'Upper GI_Null reference',...
%                  'Lower GI_Null reference',...
%                 'Upper GI_Forced cough',...
%                 'Upper GI_Procedure starts',...
%                 'Upper GI_Procedure ends',...
%                 'Upper GI_Cough',...
%                 'Upper GI_Fundal retroflexion',...
%                 'Upper GI_Throat spray given'};

% %             
%     eventList = {'Upper GI_Null reference',...
%                  'Lower GI_Null reference',...
%                 'Upper GI_Throat spray given',...
%                 'Upper GI_Forced cough'};
%             
%      eventList = {'Upper GI_Null reference',...
%                  'Lower GI_Null reference',...
%                 'Upper GI_Throat spray given',...
%                 'Upper GI_Cough',...
%                 'Upper GI_Procedure starts',...
%                 'Upper GI_Procedure ends'};
            
     eventList = {'Upper GI_Null reference',...
                 'Lower GI_Null reference',...
                'Upper GI_Cough'};
% 
%     eventList = {'Upper GI_Null reference',...
%                 'Lower GI_Null reference',...
%                 'Upper GI_Forced cough',...
%                 'Lower GI_Procedure starts',...
%                 'Lower GI_Procedure ends',...
%                 'Lower GI_Abdominal pressure',...
%                 'Lower GI_Biopsy sampling',...
%                 'Lower GI_Caecal intubation',...
%                 'Lower GI_Catheter in',...
%                 'Lower GI_Catheter out',...
%                 'Lower GI_Position changes',...
%                 'Lower GI_Diathermy pad used',...
%                 'Lower GI_Rectal insufflation/retroflexion'};

            
    upperGI_nullreffile = 'Upper GI_Null reference';
    lowerGI_nullreffile = 'Lower GI_Null reference';
            
    nFiles = size(eventList,2);

    fileList_raw = dir(dataDir);

    filterFun = @(x) regexpi(x, '^[0-9]{8}');
    temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
    fileList_raw = fileList_raw(~cellfun(@isempty,temp));

    nFilesRaw = size(fileList_raw,1);

    targetDiams = [0.3; 0.5; 0.7; 1.0; 3.0; 5.0; 10.0];
            
    allEventNames = {};
    allEventTimes = [];
    
    eventNames_diff = [];
    fullTable = [];

    for fileIdx = 1:nFilesRaw
        currentFolder = fileList_raw(fileIdx).name;

        test = sscanf(currentFolder, '%4d%2d%2d');

        Y = test(1);
        M = test(2);
        D = test(3);

        subFileList_raw = dir([dataDir, '/', currentFolder]);

        filterFun = @(x) regexpi(x, ['^', currentFolder, '_patient[0-9]{1,2}.csv']);
        temp = cellfun(filterFun, {subFileList_raw.name}, 'UniformOutput', false); 
        subFileList_raw = subFileList_raw(~cellfun(@isempty,temp));

        currentNPatients = size(subFileList_raw,1);

        for subfileIdx = 1:currentNPatients
            currentFile = subFileList_raw(subfileIdx).name;

            test = sscanf(currentFile, '%4d%2d%2d_patient%d.csv');

            P = test(4); % TODO ensure dates match

            [data, datatimes, eventTimes, eventNames, avSampleTime, otherVars, ~, ~, data_v, ~, ~, diameters, diameters_av, data_next, opTime2_next,  bg_data_next, fg_data_next, data_v_next, bg_data_v_next, fg_data_v_next] = loadAnnotatedData(Y,M,D,P,dataDir);

            if isnan(data) % no valid events detected
                continue;
            end
            
            tDiff = diff(eventTimes);
            tDiff = [tDiff; datatimes(end) - eventTimes(end)];
            
            eventNames_next = eventNames(2:end,1);
            eventNames_next = [eventNames_next;{'End of recording'}];
            
            if isempty(eventNames_diff)
                eventNames_diff = eventNames;
                eventNames_diff_next = eventNames_next;
                tDiff_all = tDiff;
            else
                eventNames_diff = cat(1, eventNames_diff, eventNames);
                eventNames_diff_next = cat(1, eventNames_diff_next, eventNames_next);
                tDiff_all = [tDiff_all; tDiff];
            end
            
%             allEventNames = cat(1,allEventNames,eventNames);
%             allEventTimes = cat(1,allEventTimes,eventTimes);

            tempIdxes = (1:size(eventTimes,1))';
            tempIdxes = kron(tempIdxes, ones(size(targetDiams,1),1));
            eventNames_rep = eventNames(tempIdxes,1);
            eventTimes_rep = eventTimes(tempIdxes,1);
            otherVars_rep = otherVars(ones(size(eventTimes_rep,1),1),:);
            
            eventMargin = 6;
            
            tempTimes = -60*eventMargin:60*eventMargin;
            startTime = min(tempTimes);
            endTime = max(tempTimes);
            timeStep = tempTimes(2) - tempTimes(1);
            tempTimes = seconds(repmat(tempTimes,size(eventTimes,1),1)) + repmat(eventTimes,1,size(tempTimes,2));
            
            [A,B] = ismember(tempTimes,datatimes);
            
            tempData_all = nan(size(tempTimes,1)*size(targetDiams,1),size(tempTimes,2));
            for k=1:size(targetDiams,1)
                
                idx = find(diameters == targetDiams(k));
                
                if (~isempty(idx))
                    tempData = nan(size(tempTimes));
                    data_temp = data(idx,:)';
                    tempData(A) = data_temp(B(A));
                    
                    tempData_all(k:size(targetDiams,1):end,:) = tempData;
                end
            end
            
            tempDiams_all = repmat(targetDiams, size(eventTimes,1),1);
            
            fullTabTemp = otherVars_rep;
            fullTabTemp.eventNames = table2cell(eventNames_rep);
            fullTabTemp.eventTimes = eventTimes_rep;
            fullTabTemp.diameters = tempDiams_all;
            temp = cat(2,fullTabTemp,array2table(tempData_all));
            
            if isempty(fullTable)
                fullTable = temp;
            else
                fullTable = cat(1,fullTable,temp);
            end
            
        end
    end
    
    % Compute gap between events for calibrating integration window
    for kk = 1:size(eventList,2)
    	file = eventList{kk};
        label = file;
        eventName = label(10:end);
        
        currentEventsValid = strcmpi(table2cell(eventNames_diff),eventName);
        
        meanTdiff(kk) = mean(tDiff_all(currentEventsValid));
        minTdiff(kk) = min(tDiff_all(currentEventsValid));
        maxTdiff(kk) = max(tDiff_all(currentEventsValid));
        
        if kk==1
            eventNames_diff_plot = eventNames_diff(currentEventsValid,:);
            tDiffPlot = tDiff_all(currentEventsValid);
        else
            eventNames_diff_plot  = [eventNames_diff_plot; eventNames_diff(currentEventsValid,:)];
            tDiffPlot = [tDiffPlot; tDiff_all(currentEventsValid)];
        end
        a = 1;
    end
  
    dataTemp = seconds(tDiffPlot);
    catsTemp = categorical(table2cell(eventNames_diff_plot));
    violinplot(dataTemp.', catsTemp.');
    xlabel('Event name');
    ylabel('Time until next event');
    title('Event windows');
    xtickangle(60);
    a = 1;
end

%% Table for soring vars
varNames = {'label', 'Nevents', 'Npatients', 'mean_n', 'std_n', 'median_n', 'lq_n', 'uq_n',  'mean_v', 'std_v', 'median_v', 'lq_v', 'uq_v', 'mu_mean_n', 'mu_std_n', 'sig_a_n', 'sig_b_n', 'mu_mean_v', 'mu_std_v', 'sig_a_v', 'sig_b_v'};

for k=1:size(varNames,2)
    if k == 1
        varTypes{k} = 'string';
    else
        varTypes{k} = 'double';
    end
end

startFileIdx = 1;
resultsTable = table('Size',[nFiles, size(varNames,2)],'VariableTypes', varTypes, 'VariableNames', varNames);


for fileIdx = 1:nFiles
    if fileIdx > 2 && fileIdx < startFileIdx
        continue;
    end
    
    if startFileIdx > 1
        if fileIdx == startFileIdx
            load(fullfile(folder,resultsFileName));
        end
    end
    
    if useCSVs
        file = fileList{fileIdx};
        label = file(1:end-4);
    
        isLowerGI = ~isempty(regexp(file, 'LowerGI_*.'));
        isUpperGI = ~isempty(regexp(file, 'UpperGI_*.'));
    else
        file = eventList{fileIdx};
        label = file;
        eventName = label(10:end);
        label = strrep(label,'/',' ');
        label = strrep(label,'\',' ');
        
        isLowerGI = ~isempty(regexp(file, 'Lower GI_*.'));
        isUpperGI = ~isempty(regexp(file, 'Upper GI_*.'));
    end
    
    if limitSize
        if lt
            label = [label, '_lt', num2str(sizeLim)];
            resultsFileName = ['resultsTable_lt', num2str(sizeLim), '.mat'];
        else
            label = [label, '_gt', num2str(sizeLim)];
            resultsFileName = ['resultsTable_gt', num2str(sizeLim), '.mat'];
        end
    else
        resultsFileName = 'resultsTable.mat';
    end
    
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
    
    
    if useCSVs
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
            [isSmoker, isSmokerCats] = processVars(otherVars.Smoker, {'yes','no', 'unknown'}, {'.*yes.*', '.*[1-9]+[0-9]*.*'}, {'.*no.*', '.*[0]+.*'}, {'.*none.*','.*N/A.*',''});

            % Mask
            usesPatientMask_raw = ismember(patientNos,[59,60,61, 64,65]); %FIX should pull this
            usesPatientMask = categorical(usesPatientMask_raw+1,1:2,{'no', 'yes'});

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
                [hiatusHernia, hiatusHerniaCats] = processVars(otherVars.HiatusHernia, {'no', 'yes', 'unknown'}, {'.*none.*', '.*[0].*', '.*no.*'}, {'.*[1-9]+[0-9]*.*', '.*yes.*', '.*massive.*'}, {'.*unknown.*','.*N/A.*',''});
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
    else
        %% Up to here - loading data
        %Sort by event first
        
        currentValidRows = strcmpi(fullTable.eventNames, eventName);
        
        if isUpperGI
            currentValidRows = currentValidRows & fullTable.Procedure == 'upper GI';
        elseif isLowerGI
            currentValidRows = currentValidRows & fullTable.Procedure == 'lower GI';
        end
        T = fullTable(currentValidRows,:);

        [tempTimes, tempIdxes] = unique(T.eventTimes);
        
        dataStartCol = find(strcmpi(T.Properties.VariableNames,'tempData_all1'));
        otherVars = T(sort(tempIdxes, 'ascend'),1:(dataStartCol-2));
        
        eventTimes = tempTimes;
        indices = 1:size(eventTimes,1);
        indices = indices';
        diameters_full = T.diameters;
        diameters = diameters_full(T.eventTimes == eventTimes(1));
        
        patientNos = otherVars.StudyNumber;

        sex = otherVars.Sex;
        sedation = otherVars.Sedation;
        analTone = otherVars.AnalTone;
        useOfCO2orWater = otherVars.UseOfCO2orWater;
        isSmoker = otherVars.Smoker;
        usesPatientMask = otherVars.PatientMask;
        roomType = otherVars.RoomType;
        ugiRoute = otherVars.UGIroute;
        diverticularDisease = otherVars.DiverticularDisease;
        looping = otherVars.Looping;
        discomfort = otherVars.Discomfort;
        hiatusHernia = otherVars.HiatusHernia;
        suctioning = otherVars.Suctioning;
        procedureType = otherVars.ProcedureType;

        % Add air sentry and toher vars
    end

    %% Correct for effect of tube
    useTubeCorrection = true;

    if useTubeCorrection
        
        if useCSVs
            tubeCorrection_tab = readtable(['C:\Users\', username, '\OneDrive - The University of Nottingham\SAVE\TubeCalibration\TubeBendCorrection.csv']);
            %tubeCorrection_tab = readtable('/home/george/Desktop/TubeBendCorrection.csv');
            tubeCorrection = table2array(tubeCorrection_tab);

            for k=1:size(diameters,1)
                correctionIdx = find(diameters(k) == tubeCorrection(:,1));

                correctionVal(k) = tubeCorrection(correctionIdx,2);
            end
        else
            % Direct method automatically applies tube correction
            correctionVal = zeros(1,size(diameters,1));
        end
    else
        correctionVal = zeros(1,size(diameters,1));
    end

    %%

    maxDiameter = 25; % Should load this from datasheet?
    diameters = [diameters; maxDiameter];

    time = startTime:timeStep:endTime;
    nTimes = (endTime - startTime)/timeStep + 1;

    if useCSVs
        validRows = table2array(T(:,4)) == eventTimes(1); % FIX should loop for other indices and check
    else
        validRows = (T.eventTimes) == eventTimes(1); % FIX should loop for other indices and check
    end
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

        if useCSVs
            validRows = table2array(T(:,4)) == eventTimes(currentIdx); 
        else
            validRows = (T.eventTimes) == eventTimes(currentIdx);
        end

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
        
        if useCSVs
            correctForAT = true; %AeroTrak time is the start of the sample. Jane's CSVs do not do this.
        else
            correctForAT = false;
        end
        if correctForAT
            data = circshift(data, avSampleTime);
            data(1:avSampleTime) = NaN;
            data_v = circshift(data_v, avSampleTime);
            data_v(1:avSampleTime) = NaN;
        end
        
        % Densities so that a probability density approach can be used
        data_v_density = data_v ./ repmat(currentLogBinSizes,1, size(data,2)); %Try using log binsizes
        data_density = data ./ repmat(currentLogBinSizes,1, size(data,2));
        
        [bg_current, fg_current] = splitBGFG(data, avSampleTime, tempValid, 'causal', false);

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
        if ~isempty(regexpi(label, '.*Forced cough.*'))
            rawdiffwindowSize_after = 21; %Small to avoid throat spray
        elseif ~isempty(regexpi(label, '.*Deep breath.*'))
            rawdiffwindowSize_after = 15; % Much shorter due to quick succession of measurements
        elseif ~isempty(regexpi(label, '.*Speaking.*'))
            rawdiffwindowSize_after = 15; % Much shorter due to quick succession of measurements
        elseif ~isempty(regexpi(label, '.*Insufflation.*'))
            rawdiffwindowSize_after = 15; % Much shorter due to quick succession of measurements
        elseif ~isempty(regexpi(label, '.*Fundal.*'))
            rawdiffwindowSize_after = 40; % Much shorter due to quick succession of measurements
        else
            rawdiffwindowSize_after = 100;
        end
        
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
                

            %else
            %    A_diff(k) = A_est;
            %    A_v_diff(k) = A_v_est;
            %end
            
            if limitSize
                lt5 = currentDiameters < sizeLim;
                lt5 = lt5(1:end-1); % To get rid of upper bin bound
                A_t_lt5 = sum(currentData_raw(lt5));
                A_t_gt5 = sum(currentData_raw(~lt5));

                if lt
                    A_diff(k) = A_t_lt5;
                else
                    A_diff(k) = A_t_gt5;
                end
            else
                A_diff(k) = A_t;
            end
            A_v_diff(k) = A_t_v;
     
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
                reject = A_diff < -4*noiseStd; % Reject points that are probably errors
                reject_v = A_v_diff < -4*noiseStd_v; % Reject points that are probably errors
                
                label_temp = label;
                
                if rejectMasks
                    if ~strcmpi(label_temp, 'Upper GI_Deep breaths') && ~strcmpi(label_temp, 'Upper GI_Forced cough') && ~strcmpi(label_temp, 'Upper GI_Speaking')
                        reject = reject | (usesPatientMask ~= 'no')';
                        reject_v = reject_v | (usesPatientMask ~= 'no')';
                    end
                    label = [label, '_rejectMasks'];
                end
                
                if rejectERCP
                    reject = reject | (procedureType == 'ERCP')' | (procedureType == 'EUS')';
                    reject_v = reject_v | (procedureType == 'ERCP')' | (procedureType == 'EUS')';
                end
                
                if rejectLF
                    reject = reject | (roomType == 'laminar flow')';
                    reject_v = reject_v | (roomType == 'laminar flow')';
                end
                
                if LFonly
                    reject = reject | (roomType ~= 'laminar flow')';
                    reject_v = reject_v | (roomType ~= 'laminar flow')';
                    
                    noToKeep = nnz(~reject)
                end
                
                if cytospongeOnly
                    reject = reject | (procedureType ~= 'cytosponge')';
                    reject_v = reject_v | (procedureType ~= 'cytosponge')';
                end
                
                if rejectCytosponge
                    reject = reject | (procedureType == 'cytosponge')';
                    reject_v = reject_v | (procedureType == 'cytosponge')';
                end

                if rejectTheatre
                    if ~strcmpi(label_temp, 'Upper GI_Deep breaths') && ~strcmpi(label_temp, 'Upper GI_Forced cough') && ~strcmpi(label_temp, 'Upper GI_Speaking')
                        reject = reject | (roomType == 'theatre')';
                        reject_v = reject_v | (roomType == 'theatre')';
                    end
                end
                
                if rejectNasalAbandoned
                    if ~strcmpi(label_temp, 'Upper GI_Deep breaths') && ~strcmpi(label_temp, 'Upper GI_Forced cough') && ~strcmpi(label_temp, 'Upper GI_Speaking')
                        reject = reject | (ugiRoute == 'nasal abandoned')';
                        reject_v = reject_v | (ugiRoute == 'nasal abandoned')';
                    end
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
                options_sa = optimoptions('simulannealbnd','MaxFunctionEvaluations',1e5, 'FunctionTolerance',1e-12, 'Display', 'off');
                %bestVal = fmincon(objFun,x0,[],[],[],[], LB , UB,[],options_fmincon);
                bestVal = simulannealbnd(objFun,x0, LB , UB, options_sa);
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
    resultsTable.sig_raw(fileIdx) = {sigma_diff(~reject)};
    resultsTable.sig_v_raw(fileIdx) = {sigma_v_diff(~reject_v)};
    
    resultsTable.age(fileIdx) = {otherVars.Age(~reject)};
    resultsTable.sex(fileIdx) = {otherVars.Sex(~reject)};
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
    resultsTable.procedureType(fileIdx) = {otherVars.ProcedureType(~reject)};
    
    resultsTable.age_v(fileIdx) = {otherVars.Age(~reject_v)};
    resultsTable.sex_v(fileIdx) = {otherVars.Sex(~reject_v)};
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
    resultsTable.procedureType_v(fileIdx) = {otherVars.ProcedureType(~reject_v)};
    
    resultsTable.samples{fileIdx} = []; %Delete
    
    % Compute pvals
    for tempIdx = 1:fileIdx-1
        event1 = resultsTable.n_raw{tempIdx};
        event1 = event1(:);
        event1_name = resultsTable.label(tempIdx);
        event2 = resultsTable.n_raw{fileIdx};
        event2 = event2(:);
        event2_name = resultsTable.label(fileIdx);
        
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
        

        
%         %REMOVE
%         if tempIdx == 19 || fileIdx == 19 %forcedCough
%             computeEventPvals = true;
%         else
%             computeEventPvals = false;
%         end
        
        if computeEventPvals && (~isempty(regexpi(event1_name, 'Upper GI_Null reference.*')) || ~isempty(regexpi(event1_name, 'Upper GI_Forced cough.*')))
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [pMu, pSig, meanRatio, ratioLowCI, ratioUpperCI, mean1_out, mean2_out, samples1out, samples2out] = computeSignificance(event1, event2, noiseMean, noiseStd, mean1, std1, mean2, std2);
            resultsTable.mean_n(tempIdx) = mean1_out;
            resultsTable.mean_n(fileIdx) = mean2_out;
        else
            pMu = 0.5;
            pSig = 0.5;
            
            meanRatio = 1;
            ratioLowCI = 1;
            ratioUpperCI = 1;
        end
            
        pMu = min([pMu, 1-pMu]);
        pSig = min([pSig, 1-pSig]);
     
        pMuTable(fileIdx, tempIdx) = pMu;
        pSigTable(fileIdx, tempIdx) = pSig;
        
        meanRatioTable(fileIdx, tempIdx) = meanRatio;
        ratioLowCITable(fileIdx, tempIdx) = ratioLowCI;
        ratioUpperCITable(fileIdx, tempIdx) = ratioUpperCI;
        
        
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
        

        if computeEventPvals_v
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [pMu_v, pSig_v, meanRatio, ratioLowCI, ratioUpperCI, mean1_out, mean2_out] = computeSignificance(event1, event2, noiseMean_v, noiseStd_v, mean1, std1, mean2, std2, 'muMinIn', -34);
            resultsTable.mean_v(tempIdx) = mean1_out;
            resultsTable.mean_v(fileIdx) = mean2_out;
        else
            pMu_v = 0.5;
            pSig_v = 0.5;
        end
            
        pMu_v = min([pMu_v, 1-pMu_v]);
        pSig_v = min([pSig_v, 1-pSig_v]);
     
        pMuTable_v(fileIdx, tempIdx) = pMu_v;
        pSigTable_v(fileIdx, tempIdx) = pSig_v;
        
        
        meanRatioTable_v(fileIdx, tempIdx) = meanRatio;
        ratioLowCITable_v(fileIdx, tempIdx) = ratioLowCI;
        ratioUpperCITable_v(fileIdx, tempIdx) = ratioUpperCI;
        
        % Now mean diameter
        event1 = resultsTable.mu_raw{tempIdx};
        event1 = event1(:);
        %event1 = log(event1); % Sizes already in log so no required.
        event2 = resultsTable.mu_raw{fileIdx};
        event2 = event2(:);
        %event2 = log(event2);
               

        if computeEventPvals_mu
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [~,pMu_mu] = ttest2(event1, event2);
        else
            pMu_mu = 0.5;
        end
            
        pMu_mu = min([pMu_mu, 1-pMu_mu]);
     
        pMuTable_mu(fileIdx, tempIdx) = pMu_mu;
        
        % Now spread of sizes
        event1 = resultsTable.sig_raw{tempIdx};
        event1 = event1(:);
        event2 = resultsTable.sig_raw{fileIdx};
        event2 = event2(:);
               
    
        if computeEventPvals_sig
            disp(['Computing significance..', num2str(tempIdx), '/', num2str(fileIdx)]);
            [~,pMu_sig] = kstest2(event1, event2); %Possibly not valid for gamma distributed vars!
        else
            pMu_sig = 0.5;
        end
            
        pMu_sig = min([pMu_sig, 1-pMu_sig]);
     
        pMuTable_sig(fileIdx, tempIdx) = pMu_sig;
        
    end
    
    data_box = [];
    data_v_box = [];
    mu_box = [];
    mu_v_box = [];
    sig_box = [];
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
             
        currentRaw_sig_n = resultsTable.sig_raw{dataSetIdx};
        currentRaw_sig_n = currentRaw_sig_n(:);
        
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
            sig_box = [sig_box; currentRaw_sig_n;];
        else
            mu_box = [mu_box; zeros(size(currentRaw_mu_n));];
            mu_v_box = [mu_v_box; zeros(size(currentRaw_mu_v));];
            sig_box = [sig_box; zeros(size(currentRaw_sig_n));];
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
        
        boxFig_sig_n = figure('units','normalized','outerposition',[0 0 1 1]);
        violinFig_sig_n = figure('units','normalized','outerposition',[0 0 1 1]);
        
        
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
        dataCompFigs(5).plotPvals = true;
        dataCompFigs(5).pValVarMu = 'pMuTable_mu';
        dataCompFigs(5).pValVarSig = 'pMuTable_mu';
        dataCompFigs(5).plotNonLogYticks = true;
        dataCompFigs(5).rawData = 'mu_raw';
        dataCompFigs(5).meanVar = 'mu_mean_n';
        dataCompFigs(5).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(5).saveName = ['diameter_number_', dataCompFigs(5).plotType];
        
        dataCompFigs(6).fig = violinFig_mu_n;
        dataCompFigs(6).plotType = 'violinplot';
        dataCompFigs(6).varNames = 'mu_box';
        dataCompFigs(6).catNames = 'cats_box';
        dataCompFigs(6).plotPvals = true;
        dataCompFigs(6).pValVarMu = 'pMuTable_mu';
        dataCompFigs(6).pValVarSig = 'pMuTable_mu';
        dataCompFigs(6).plotNonLogYticks = true;
        dataCompFigs(6).rawData = 'mu_raw';
        dataCompFigs(6).meanVar = 'mu_mean_n';
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
        dataCompFigs(7).meanVar = 'mu_mean_v';
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
        dataCompFigs(8).meanVar = 'mu_mean_v';
        dataCompFigs(8).ylabel = 'Mean particle diameter (\mum)';
        dataCompFigs(8).saveName = ['diameter_volume_', dataCompFigs(8).plotType];
        
        dataCompFigs(9).fig = boxFig_sig_n;
        dataCompFigs(9).plotType = 'boxplot';
        dataCompFigs(9).varNames = 'sig_box';
        dataCompFigs(9).catNames = 'cats_box';
        dataCompFigs(9).plotPvals = true;
        dataCompFigs(9).pValVarMu = 'pMuTable_sig';
        dataCompFigs(9).pValVarSig = 'pMuTable_sig';
        dataCompFigs(9).plotNonLogYticks = true;
        dataCompFigs(9).rawData = 'sig_raw';
        dataCompFigs(9).meanVar = '';
        dataCompFigs(9).ylabel = 'particle diameter spread (\mum)';
        dataCompFigs(9).saveName = ['diameterspread_number_', dataCompFigs(5).plotType];
        
        dataCompFigs(10).fig = violinFig_sig_n;
        dataCompFigs(10).plotType = 'violinplot';
        dataCompFigs(10).varNames = 'sig_box';
        dataCompFigs(10).catNames = 'cats_box';
        dataCompFigs(10).plotPvals = true;
        dataCompFigs(10).pValVarMu = 'pMuTable_sig';
        dataCompFigs(10).pValVarSig = 'pMuTable_sig';
        dataCompFigs(10).plotNonLogYticks = true;
        dataCompFigs(10).rawData = 'sig_raw';
        dataCompFigs(10).meanVar = '';
        dataCompFigs(10).ylabel = 'particle diameter spread (\mum)';
        dataCompFigs(10).saveName = ['diameterspread_number_', dataCompFigs(6).plotType];
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
                            %computeAndPlotPvals(tempD1,tempD2,[],[],tempIdx, fileIdx2,'pMeanIn', pValMu, 'pSigIn', pValSig);
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
    
    
    if analyseVariables %% Now plot variables
        if (fileIdx == 1)
            sedationFig_n = figure;
            sexFig_n = figure;
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
            procedureTypeFig_n = figure;
            ageFig_n = figure;
            bmiFig_n = figure;

            sedationFig_v = figure;
            sexFig_v = figure;
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
            procedureTypeFig_v = figure;
            ageFig_v = figure;
            bmiFig_v = figure;

            sedationFig_mu = figure;
            sexFig_mu = figure;
            analToneFig_mu = figure;
            co2vWaterFig_mu = figure;
            smokerFig_mu = figure;
            maskFig_mu = figure;
            roomTypeFig_mu = figure;
            ugiRouteFig_mu = figure;
            divertDiseaseFig_mu = figure;
            loopingFig_mu = figure;
            discomfortFig_mu = figure;
            hiatusFig_mu = figure;
            suctioningFig_mu = figure;
            procedureTypeFig_mu = figure;
            ageFig_mu = figure;
            bmiFig_mu = figure;

            sedationFig_sig = figure;
            sexFig_sig = figure;
            analToneFig_sig = figure;
            co2vWaterFig_sig = figure;
            smokerFig_sig = figure;
            maskFig_sig = figure;
            roomTypeFig_sig = figure;
            ugiRouteFig_sig = figure;
            divertDiseaseFig_sig = figure;
            loopingFig_sig = figure;
            discomfortFig_sig = figure;
            hiatusFig_sig = figure;
            suctioningFig_sig = figure;
            procedureTypeFig_sig = figure;
            ageFig_sig = figure;
            bmiFig_sig = figure;

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

            %repeat for size
            varFigs(29).fig = sedationFig_mu;
            varFigs(29).varname = 'sedation';
            varFigs(29).rawvar = 'mu_raw';
            varFigs(29).type = 'discrete';
            varFigs(30).fig = analToneFig_mu;
            varFigs(30).varname = 'analTone';
            varFigs(30).rawvar = 'mu_raw';
            varFigs(30).type = 'discrete';
            varFigs(31).fig = co2vWaterFig_mu;
            varFigs(31).varname = 'useOfCO2OrWater';
            varFigs(31).rawvar = 'mu_raw';
            varFigs(31).type = 'discrete';
            varFigs(32).fig = smokerFig_mu;
            varFigs(32).varname = 'isSmoker';
            varFigs(32).rawvar = 'mu_raw';
            varFigs(32).type = 'discrete';
            varFigs(33).fig = maskFig_mu;
            varFigs(33).varname = 'usesPatientMask';
            varFigs(33).rawvar = 'mu_raw';
            varFigs(33).type = 'discrete';
            varFigs(34).fig = roomTypeFig_mu;
            varFigs(34).varname = 'roomType';
            varFigs(34).rawvar = 'mu_raw';
            varFigs(34).type = 'discrete';
            varFigs(35).fig = ugiRouteFig_mu;
            varFigs(35).varname = 'ugiRoute';
            varFigs(35).rawvar = 'mu_raw';
            varFigs(35).type = 'discrete';
            varFigs(36).fig = divertDiseaseFig_mu;
            varFigs(36).varname = 'diverticularDisease';
            varFigs(36).rawvar = 'mu_raw';
            varFigs(36).type = 'discrete';
            varFigs(37).fig = loopingFig_mu;
            varFigs(37).varname = 'looping';
            varFigs(37).rawvar = 'mu_raw';
            varFigs(37).type = 'discrete';
            varFigs(38).fig = discomfortFig_mu;
            varFigs(38).varname = 'discomfort';
            varFigs(38).rawvar = 'mu_raw';
            varFigs(38).type = 'discrete';
            varFigs(39).fig = hiatusFig_mu;
            varFigs(39).varname = 'hiatusHernia';
            varFigs(39).rawvar = 'mu_raw';
            varFigs(39).type = 'discrete';
            varFigs(40).fig = suctioningFig_mu;
            varFigs(40).varname = 'suctioning';
            varFigs(40).rawvar = 'mu_raw';
            varFigs(40).type = 'discrete';
            varFigs(41).fig = ageFig_mu;
            varFigs(41).varname = 'age';
            varFigs(41).rawvar = 'mu_raw';
            varFigs(41).type = 'continuous';
            varFigs(42).fig = bmiFig_mu;
            varFigs(42).varname = 'bmi';
            varFigs(42).rawvar = 'mu_raw';
            varFigs(42).type = 'continuous';

            %repeat for size spread
            varFigs(43).fig = sedationFig_sig;
            varFigs(43).varname = 'sedation';
            varFigs(43).rawvar = 'sig_raw';
            varFigs(43).type = 'discrete';
            varFigs(44).fig = analToneFig_sig;
            varFigs(44).varname = 'analTone';
            varFigs(44).rawvar = 'sig_raw';
            varFigs(44).type = 'discrete';
            varFigs(45).fig = co2vWaterFig_sig;
            varFigs(45).varname = 'useOfCO2OrWater';
            varFigs(45).rawvar = 'sig_raw';
            varFigs(45).type = 'discrete';
            varFigs(46).fig = smokerFig_sig;
            varFigs(46).varname = 'isSmoker';
            varFigs(46).rawvar = 'sig_raw';
            varFigs(46).type = 'discrete';
            varFigs(47).fig = maskFig_sig;
            varFigs(47).varname = 'usesPatientMask';
            varFigs(47).rawvar = 'sig_raw';
            varFigs(47).type = 'discrete';
            varFigs(48).fig = roomTypeFig_sig;
            varFigs(48).varname = 'roomType';
            varFigs(48).rawvar = 'sig_raw';
            varFigs(48).type = 'discrete';
            varFigs(49).fig = ugiRouteFig_sig;
            varFigs(49).varname = 'ugiRoute';
            varFigs(49).rawvar = 'sig_raw';
            varFigs(49).type = 'discrete';
            varFigs(50).fig = divertDiseaseFig_sig;
            varFigs(50).varname = 'diverticularDisease';
            varFigs(50).rawvar = 'sig_raw';
            varFigs(50).type = 'discrete';
            varFigs(51).fig = loopingFig_sig;
            varFigs(51).varname = 'looping';
            varFigs(51).rawvar = 'sig_raw';
            varFigs(51).type = 'discrete';
            varFigs(52).fig = discomfortFig_sig;
            varFigs(52).varname = 'discomfort';
            varFigs(52).rawvar = 'sig_raw';
            varFigs(52).type = 'discrete';
            varFigs(53).fig = hiatusFig_sig;
            varFigs(53).varname = 'hiatusHernia';
            varFigs(53).rawvar = 'sig_raw';
            varFigs(53).type = 'discrete';
            varFigs(54).fig = suctioningFig_sig;
            varFigs(54).varname = 'suctioning';
            varFigs(54).rawvar = 'sig_raw';
            varFigs(54).type = 'discrete';
            varFigs(55).fig = ageFig_sig;
            varFigs(55).varname = 'age';
            varFigs(55).rawvar = 'sig_raw';
            varFigs(55).type = 'continuous';
            varFigs(56).fig = bmiFig_sig;
            varFigs(56).varname = 'bmi';
            varFigs(56).rawvar = 'sig_raw';
            varFigs(56).type = 'continuous';

            % Sex - variable I forgot!
            varFigs(57).fig = sexFig_n;
            varFigs(57).varname = 'sex';
            varFigs(57).rawvar = 'n_raw';
            varFigs(57).type = 'discrete';

            varFigs(58).fig = sexFig_v;
            varFigs(58).varname = 'sex_v';
            varFigs(58).rawvar = 'v_raw';
            varFigs(58).type = 'discrete';

            varFigs(59).fig = sexFig_mu;
            varFigs(59).varname = 'sex';
            varFigs(59).rawvar = 'mu_raw';
            varFigs(59).type = 'discrete';

            varFigs(60).fig = sexFig_sig;
            varFigs(60).varname = 'sex';
            varFigs(60).rawvar = 'sig_raw';
            varFigs(60).type = 'discrete';
            
            %Procedure type - variable I forgot!
            varFigs(61).fig = procedureTypeFig_n;
            varFigs(61).varname = 'procedureType';
            varFigs(61).rawvar = 'n_raw';
            varFigs(61).type = 'discrete';

            varFigs(62).fig = procedureTypeFig_v;
            varFigs(62).varname = 'procedureType_v';
            varFigs(62).rawvar = 'v_raw';
            varFigs(62).type = 'discrete';

            varFigs(63).fig = procedureTypeFig_mu;
            varFigs(63).varname = 'procedureType';
            varFigs(63).rawvar = 'mu_raw';
            varFigs(63).type = 'discrete';

            varFigs(64).fig = procedureTypeFig_sig;
            varFigs(64).varname = 'procedureType';
            varFigs(64).rawvar = 'sig_raw';
            varFigs(64).type = 'discrete';
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
                    if strcmpi(varFigs(k).rawvar, 'n_raw') || strcmpi(varFigs(k).rawvar, 'mu_raw') || strcmpi(varFigs(k).rawvar, 'sig_raw')
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

                         if fileIdx > 2 && strcmpi(varFigs(k).varname, 'procedureType') %% So that it only does cytosponge
                            if computeVarPvals
                                if strcmpi(varFigs(k).rawvar, 'n_raw')
                                    if LFonly
                                        [~, pValMu] = ttest2(log(d1),log(d2));
                                        computeAndPlotPvals(d1,d1,[],[],cat1Idx, cat2Idx,'pMeanIn', pValMu, 'pSigIn', pValMu);
                                    else
                                        computeAndPlotPvals(d1,d2,noiseMean,noiseStd,cat1Idx,cat2Idx, 'pThresh', 0.5);
                                    end
                                elseif strcmpi(varFigs(k).rawvar, 'v_raw')
                                    %computeAndPlotPvals(d1,d2,noiseMean,noiseStd,cat1Idx,cat2Idx, 'pThresh', 0.5, 'muMinIn', -34);
                                elseif strcmpi(varFigs(k).rawvar, 'mu_raw')
                                    [~, pValMu] = ttest2(d1,d2);
                                    computeAndPlotPvals(d1,d1,[],[],cat1Idx, cat2Idx,'pMeanIn', pValMu, 'pSigIn', pValMu, 'pThresh', 0.5);
                                elseif strcmpi(varFigs(k).rawvar, 'sig_raw')
                                    [~, pValSig] = kstest2(d1,d2);
                                    computeAndPlotPvals(d1,d1,[],[],cat1Idx, cat2Idx,'pMeanIn', pValSig, 'pSigIn', pValSig);
                                end
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
            if strcmpi(varFigs(k).rawvar, 'n_raw')  || strcmpi(varFigs(k).rawvar, 'mu_raw') || strcmpi(varFigs(k).rawvar, 'sig_raw')
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
        elseif strcmpi(varFigs(k).rawvar, 'mu_raw')
            ylabel('Mean particle diameter (log um)');
        elseif strcmpi(varFigs(k).rawvar, 'sig_raw')
            ylabel('Mean particle spread (log um)');
        end
        xlabel(varFigs(k).varname);
        title(resultsTable.label(fileIdx), 'Interpreter', 'none');
        
        if saveFigs
            if strcmpi(varFigs(k).rawvar, 'mu_raw')
                saveFigName = [label, '_', varFigs(k).varname, '_size'];
            elseif strcmpi(varFigs(k).rawvar, 'sig_raw')
                saveFigName = [label, '_', varFigs(k).varname, '_spread'];
            else
                saveFigName = [label, '_', varFigs(k).varname];
            end
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
    end
    
    %%
    
    if (fileIdx == 1)
        dTree_fig = figure;
        forest_fig = figure('units','normalized','outerposition',[0 0 1 1]);
    end
        
    if (fileIdx > 2)
    
        % Make a table to fit
        n_raw_tab = resultsTable.n_raw{fileIdx}';
        v_raw_tab = resultsTable.v_raw{fileIdx}';
        mu_raw_tab = resultsTable.mu_raw{fileIdx}';
        sig_raw_tab = resultsTable.sig_raw{fileIdx}';
        mu_v_raw_tab = resultsTable.mu_v_raw{fileIdx}';
        age_tab = resultsTable.age{fileIdx};
        sex_tab = resultsTable.sex{fileIdx};
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

        tempTab = table(age_tab, sex_tab, bmi_tab, sedation_tab, analTone_tab, co2water_tab, isSmoker_tab, patientMask_tab, roomType_tab, ugiRoute_tab, divertDisease_tab, looping_tab, discomfort_tab, hiatusHernia_tab, suctioning_tab);

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
        
        tempForest = TreeBagger(1000, tempTab, n_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on', 'PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

        %view(tempForest.Trees{1},'Mode','graph')
        imp_perm = tempForest.OOBPermutedPredictorDeltaError;
        imp = tempForest.DeltaCriterionDecisionSplit;
        imp_mat = tempForest.SurrogateAssociation;
        resultsTable.predictorImportance(fileIdx) = {imp};
        resultsTable.predictorImportance_perm(fileIdx) = {imp_perm};
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,3,1);
        bar(imp);
        title('Variable importance');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,3,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        subplot(1,3,3);
        imagesc(imp_mat);
        axis equal;
        axis tight
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';
        h.YTickLabel = tempForest.PredictorNames;
        xticks(1:size(tempForest.PredictorNames,2));
        yticks(1:size(tempForest.PredictorNames,2));
        
        if saveFigs
            saveFigName = [label, '_random_forest'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
        
        %% Now for volumes
        age_tab_v = resultsTable.age_v{fileIdx};
        sex_tab_v = resultsTable.sex_v{fileIdx};
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

        tempTab_v = table(age_tab_v, sex_tab_v, bmi_tab_v, sedation_tab_v, analTone_tab_v, co2water_tab_v, isSmoker_tab_v, patientMask_tab_v, roomType_tab_v, ugiRoute_tab_v, divertDisease_tab_v, looping_tab_v, discomfort_tab_v, hiatusHernia_tab_v, suctioning_tab_v);

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
        
        tempForest = TreeBagger(1000, tempTab_v, v_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on', 'PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

        %view(tempForest.Trees{1},'Mode','graph')
        imp_perm = tempForest.OOBPermutedPredictorDeltaError;
        imp = tempForest.DeltaCriterionDecisionSplit;
        imp_mat = tempForest.SurrogateAssociation;
        resultsTable.predictorImportance_v(fileIdx) = {imp};
        resultsTable.predictorImportance_perm_v(fileIdx) = {imp_perm};
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,3,1);
        bar(imp);
        title('Variable importance vol.');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,3,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        subplot(1,3,3);
        imagesc(imp_mat);
        axis equal;
        axis tight
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';
        h.YTickLabel = tempForest.PredictorNames;
        xticks(1:size(tempForest.PredictorNames,2));
        yticks(1:size(tempForest.PredictorNames,2));
        
        
        if saveFigs
            saveFigName = [label, '_random_forest_v'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
        
        %% Now for sizes
        tempTree = fitrtree(tempTab, mu_raw_tab, 'CategoricalPredictors',[3:size(table,2)],'Surrogate','on');
        before = findall(groot,'Type','figure'); % Find all figures
        view(tempTree,'Mode','graph');
        after = findall(groot,'Type','figure');
        h = setdiff(after,before); % Get the figure handle of the tree viewer
        if saveFigs
            saveFigName = [label, '_regression_tree_mu'];
            saveas(h,fullfile(folder,[saveFigName, '.fig']));
            saveas(h,fullfile(folder,[saveFigName, '.png']));
        end
        close(h);
        
        tempForest = TreeBagger(1000, tempTab, mu_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on', 'PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

        %view(tempForest.Trees{1},'Mode','graph')
        imp_perm = tempForest.OOBPermutedPredictorDeltaError;
        imp = tempForest.DeltaCriterionDecisionSplit;
        imp_mat = tempForest.SurrogateAssociation;
        resultsTable.predictorImportance_mu(fileIdx) = {imp};
        resultsTable.predictorImportance_perm_mu(fileIdx) = {imp_perm};
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,3,1);
        bar(imp);
        title('Variable importance particle size.');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,3,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        subplot(1,3,3);
        imagesc(imp_mat);
        axis equal;
        axis tight
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';
        h.YTickLabel = tempForest.PredictorNames;
        xticks(1:size(tempForest.PredictorNames,2));
        yticks(1:size(tempForest.PredictorNames,2));
        
        if saveFigs
            saveFigName = [label, '_random_forest_mu'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
        
        %% Now for size spread
        tempTree = fitrtree(tempTab, sig_raw_tab, 'CategoricalPredictors',[3:size(table,2)],'Surrogate','on');
        before = findall(groot,'Type','figure'); % Find all figures
        view(tempTree,'Mode','graph');
        after = findall(groot,'Type','figure');
        h = setdiff(after,before); % Get the figure handle of the tree viewer
        if saveFigs
            saveFigName = [label, '_regression_tree_sig'];
            saveas(h,fullfile(folder,[saveFigName, '.fig']));
            saveas(h,fullfile(folder,[saveFigName, '.png']));
        end
        close(h);
        
        tempForest = TreeBagger(1000, tempTab, sig_raw_tab, 'CategoricalPredictors',[3:size(table,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on', 'PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

        %view(tempForest.Trees{1},'Mode','graph')
        imp_perm = tempForest.OOBPermutedPredictorDeltaError;
        imp = tempForest.DeltaCriterionDecisionSplit;
        imp_mat = tempForest.SurrogateAssociation;
        resultsTable.predictorImportance_sig(fileIdx) = {imp};
        resultsTable.predictorImportance_perm_sig(fileIdx) = {imp_perm};
        %imp(imp < 0) = 0;
        figure(forest_fig);
        subplot(1,3,1);
        bar(imp);
        title('Variable importance particle spread.');
        ylabel('Predictor importance estimates');
        xlabel('Predictors');
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';

        subplot(1,3,2);
        oobErrorBaggedEnsemble = oobError(tempForest);
        plot(oobErrorBaggedEnsemble)
        xlabel('Number of grown trees');
        ylabel('Out-of-bag MSE');
        
        subplot(1,3,3);
        imagesc(imp_mat);
        axis equal;
        axis tight
        h = gca;
        h.XTickLabel = tempForest.PredictorNames;
        h.XTickLabelRotation = 45;
        h.TickLabelInterpreter = 'none';
        h.YTickLabel = tempForest.PredictorNames;
        xticks(1:size(tempForest.PredictorNames,2));
        yticks(1:size(tempForest.PredictorNames,2));
        
        if saveFigs
            saveFigName = [label, '_random_forest_spread'];
            saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
            saveas(gcf,fullfile(folder,[saveFigName, '.png']));
        end
    end
    
    end

    %% Save results table
    if fileIdx >= startFileIdx
        if exist('pMuTable')
            save(fullfile(folder,resultsFileName), 'resultsTable', 'pMuTable', 'pSigTable', 'pMuTable_v', 'pSigTable_v', 'pMuTable_mu', 'pMuTable_sig', 'meanRatioTable', 'ratioLowCITable', 'ratioUpperCITable', 'meanRatioTable_v', 'ratioLowCITable_v', 'ratioUpperCITable_v');
        else
            save(fullfile(folder,resultsFileName), 'resultsTable');
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
