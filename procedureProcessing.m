%% Process individual procedures
% Author: George Gordon
% Date: 18/12/2020


% Clear any old data
clc;
clear variables;
close all;

username = getenv('username');
dataDir = ['C:\Users\', username, '\OneDrive - The University of Nottingham\SAVE\'];

fileList_raw = dir(dataDir);
    
filterFun = @(x) regexpi(x, '^[0-9]{8}');
%filterFun = @(x) regexpi(x, '^2022[0-9]{4}');
temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
fileList_raw = fileList_raw(~cellfun(@isempty,temp));

nFiles = size(fileList_raw,1);

procedureTable = [];
eventNames_all = [];
eventTimes_all = [];
interProcedureTable = [];

limitSize = true;
lt = true;
sizeLim = 5;

cytospongeOnly = false;
excludeCytosponge = false;
excludeLaminarFlow = true;
excludeTheatre  = true;
includeThroatSpray = true; % This doesn't really suppress throatspray, only removes it as the first event
excludeThroatspray = false;
suppressPosChanges = true;
suppressThroatSpray = false;
excludeMasks = true;
fgOnly = true;

%resultsFolder = sprintf('%s_results', datestr(now,'mm-dd-yyyy'));
resultsFolder = sprintf('%s_results', datestr(now,'yyyy-mm-dd'));
%folder = '01-07-2021_results';

myDiaryFile = [dataDir, resultsFolder, '/cmdOutput.txt'];

resultsFolder = [dataDir, resultsFolder];
if ~exist(resultsFolder)
    mkdir(resultsFolder);
end

diary(myDiaryFile);

label = '';
if limitSize
    if lt
        label = [label, '_lt', num2str(sizeLim)];
    else
        label = [label, '_gt', num2str(sizeLim)];
    end
end

if includeThroatSpray
    label = [label, '_ts'];
end

if fgOnly
    label = [label, '_fg'];
end

if suppressPosChanges
    label = [label, '_supPos'];
end

if suppressThroatSpray
    label = [label, '_supTS'];
end

procedureGroup = 1;
for fileIdx = 1:nFiles
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
        disp(currentFile);
        
        P = test(4); % TODO ensure dates match
        
        if ~(M==10 && D==24)% && P==5)
            %continue;
        end
        
        [data, datatimes, eventTimes, eventNames, avSampleTime, otherVars, bg_data, fg_data, data_v, bg_data_v, fg_data_v, diameters, diameters_av, data_next, opTime2_next,  bg_data_next, fg_data_next, data_v_next, bg_data_v_next, fg_data_v_next] = loadAnnotatedData(Y,M,D,P, dataDir);
        
        if istable(eventNames)
            eventNames_all = [eventNames_all; eventNames];
            eventTimes_all = [eventTimes_all; eventTimes];
        end
        
        if fgOnly
            dataSource = fg_data;
            dataSource_v = fg_data_v;
            dataSource(dataSource < 0) = 0;
            dataSource_v(dataSource_v < 0) = 0;
        else
            dataSource = data;
            dataSource_v = data_v;
        end
        
        if isnan(data) % no valid events detected
            continue;
        end
        
        if cytospongeOnly
            if otherVars.ProcedureType == 'cytosponge'
            else
                continue;
            end
        end
        
        if excludeCytosponge
            if otherVars.ProcedureType == 'cytosponge'
                continue;
            else
            end
        end
        
        if (otherVars.ProcedureType == 'ERCP') || (otherVars.ProcedureType == 'EUS')   || (otherVars.UGIroute == 'nasal abandoned')
            %Exclude masks and ERCP procedures
            continue;
        end
        
        if excludeMasks
            if (otherVars.PatientMask ~= 'no') 
                continue;
            end
        end

        if eventTimes(1) == datetime(2020,10,13)
            %Exclude ERCP
            continue;
        end

        if excludeLaminarFlow
            if (otherVars.RoomType == 'laminar flow')
               %Only use laminar flow
              continue;
            end
        end
        
        if excludeThroatspray
            if (otherVars.Sedation == 'throat spray')
               %Only use laminar flow
              continue;
            end
        end
        
        if excludeTheatre
            if (otherVars.RoomType == 'theatre')
                continue
            end
        end
        
        if any(strcmpi(table2cell(eventNames), 'forced cough'))
            a = 1;
        end
        %% Copy files
        %datestring_temp = sprintf('%0.4d%0.2d%0.2d',Y,M,D);
        %folder_src = ['C:\Users\george\OneDrive - The University of Nottingham\SAVE\',datestring_temp,'\'];
        
        %mkdir(['.\', datestring_temp]);
        %copyfile([folder_src, currentFile],['.\', datestring_temp, '\',currentFile]);
        %copyfile([folder_src, datestring_temp,'_aerotrak.xlsx'],['.\', datestring_temp, '\',datestring_temp,'_aerotrak.xlsx']);
        
        if limitSize
            if lt
                validDiams = diameters < sizeLim;
            else
                validDiams = diameters >= sizeLim;
            end
            
        else
            validDiams = true(size(diameters));
        end
        validDiams = validDiams(1:end-1); % To get rid of upper bin bound
     
        
        % Reference data
        firstEventIdx = find(datatimes >= eventTimes(1), 1);
        
        if (firstEventIdx == 1) % Not enough time before patient enters room, use another sample instead (halfway between procedure start
            procStartEventIdx = find(strcmpi(table2cell(eventNames), 'procedure starts'),1);
            procStartTime = eventTimes(procStartEventIdx);
            fakeFirstEventTime = datatimes(1) + (procStartTime - datatimes(1))/2;
            firstEventIdx = find(datatimes >= fakeFirstEventTime, 1);
        end
            
        preProcedure = dataSource(:,1:firstEventIdx-1);
        preProcedure_v = dataSource(:,1:firstEventIdx-1);
        preProcedure_av = mean(nansum(preProcedure(validDiams,:),1),2);
        [preProcedure_max, ~] = max(preProcedure,[],2);
        [preProcedure_max_sum, ~] = max(nansum(preProcedure(validDiams,:),1),[],2);
        [preProcedure_max_sum_v, ~] = max(nansum(preProcedure_v(validDiams,:),1),[],2);
        
        % Identify first event
        firstEvent = eventNames(1,1);
        firstEvent = table2cell(firstEvent);
        firstEvent = firstEvent{1};
    
        if ~isempty(regexpi(firstEvent, '.*doors locked.*')) || ~isempty(regexpi(firstEvent, '.*patient enters.*'))
            %disp(['Event ', num2str(size(procedureTable,1)+1), ': name: ', firstEvent{1}]);
            %plot(nansum(preProcedure(validDiams,:),1));
            %hold on;
            %pause(0.1);
            validPreEvent(size(procedureTable,1)+1,1) = 1;
        else
            validPreEvent(size(procedureTable,1)+1,1) = 0;
        end
%         
%         a = 1;
%         filterFun = @(x) regexpi(x, '.*prc.*');
%         temp = cellfun(filterFun, table2cell(eventNames), 'UniformOutput', false); 
%         temp2 = any(~cellfun(@isempty,temp));
%         if temp2
%             a = 1;
%         end


        
        preWindow = 0;
        postWindow = 1.4;
        
        %compFun = @(x) strcmpi(x,'procedure starts');
        %procStart_temp = find(cellfun(compFun,table2cell(eventNames)));

        compFun = @(x) regexpi(x,'.*procedure starts.*');
        procStart_temp_ = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
        procStart_temp_ = cellfun(@isempty,procStart_temp_);
        procStart_temp_ = ~procStart_temp_;
        procStart_temp = find(procStart_temp_);
        
        if isempty(procStart_temp)
            procStartIdx = 1;
            procStartIdx_raw = procStartIdx;
        else
            procStart_temp = procStart_temp(1); %somtimes theres are 2 procedure ends!
            procStartIdx = find(datatimes >= (eventTimes(procStart_temp) - minutes(preWindow)), 1);
            procStartIdx_raw = find(datatimes >= (eventTimes(procStart_temp)), 1);
        end
        
        if includeThroatSpray
            compFun = @(x) regexpi(x,'.*spray given');
            throatSpray_temp_ = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
            throatSpray_temp_ = cellfun(@isempty,throatSpray_temp_);
            throatSpray_temp_ = ~throatSpray_temp_;
            throatSpray_temp = find(throatSpray_temp_);
            if isempty(throatSpray_temp) % If there is a throatspray event
                procFirstIdx =  procStartIdx;
            else
                throatSpray_temp = throatSpray_temp(1); %somtimes theres are 2 procedure ends!
                procFirstIdx = find(datatimes >= (eventTimes(throatSpray_temp) - minutes(preWindow)), 1);
            end
        else
            procFirstIdx = procStartIdx;
        end
        
        
        %compFun = @(x) strcmpi(x,'procedure ends');
        %procEnd_temp = find(cellfun(compFun,table2cell(eventNames)));
        compFun = @(x) regexpi(x,'.*procedure ends.*');
        procEnd_temp_ = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
        procEnd_temp_ = cellfun(@isempty,procEnd_temp_);
        procEnd_temp_ = ~procEnd_temp_;
        procEnd_temp = find(procEnd_temp_);
        
        if isempty(procEnd_temp)
            procEndIdx = size(dataSource,2);
            procEndIdx_raw = procEndIdx;
        else
            procEnd_temp = procEnd_temp(1); %somtimes theres are 2 procedure ends!
            procEndIdx_raw = find(datatimes >= eventTimes(procEnd_temp), 1);
            procEndIdx = find(datatimes >= eventTimes(procEnd_temp) + minutes(postWindow), 1);
            
            if isempty(procEndIdx)
                procEndIdx = size(datatimes,1);
            end
        end
        
%         compFun = @(x) strcmpi(x,'oesophageal extubation');
%         t = find(cellfun(compFun,table2cell(eventNames)));
%         if (~isempty(t))
%             a = 1;
%         end

        % Suppress position changes
        if suppressPosChanges
            compFun = @(x) regexpi(x,'.*position changes');
            posChange_temp = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
            posChange_temp = cellfun(@isempty,posChange_temp);
            posChange_temp = ~posChange_temp;
            posChange_temp = find(posChange_temp);

            suppressionWindow = 100;
            for posChange = posChange_temp.'
                startPosChIdx = find(datatimes >= (eventTimes(posChange)), 1);
                endPosChIdx = find(datatimes > (eventTimes(posChange) + seconds(suppressionWindow)), 1) - 1;
                
                endRefIdx = startPosChIdx - 5;
                startRefIdx = find(datatimes > (eventTimes(posChange) - seconds(suppressionWindow)), 1);
                
                beforeAv = mean(dataSource(:,startRefIdx:endRefIdx),2);
                beforeAv_v = mean(dataSource_v(:,startRefIdx:endRefIdx),2);

                dataSource(:,startPosChIdx:endPosChIdx) = repmat(beforeAv,1,endPosChIdx-startPosChIdx+1);
                dataSource_v(:,startPosChIdx:endPosChIdx) = repmat(beforeAv_v,1,endPosChIdx-startPosChIdx+1);
            end
            
            if ~isempty(posChange_temp)
                posChangesSuppressed = true;
            else
                posChangesSuppressed = false;
            end
        end
        
        
        if suppressThroatSpray
            compFun = @(x) regexpi(x,'.*spray given');
            throatSpray_temp = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
            throatSpray_temp = cellfun(@isempty,throatSpray_temp);
            throatSpray_temp = ~throatSpray_temp;
            throatSpray_temp = find(throatSpray_temp);

            suppressionWindow = 100;
            for throatSpray = throatSpray_temp.'
                startThSprIdx = find(datatimes >= (eventTimes(throatSpray)), 1);
                endThSprIdx = find(datatimes > (eventTimes(throatSpray) + seconds(suppressionWindow)), 1) - 1;
                
                endRefIdx = startThSprIdx - 5;
                startRefIdx = find(datatimes > (eventTimes(throatSpray) - seconds(suppressionWindow)), 1);
                
                beforeAv = mean(dataSource(:,startRefIdx:endRefIdx),2);
                beforeAv_v = mean(dataSource_v(:,startRefIdx:endRefIdx),2);

                dataSource(:,startThSprIdx:endThSprIdx) = repmat(beforeAv,1,endThSprIdx-startThSprIdx+1);
                dataSource_v(:,startThSprIdx:endThSprIdx) = repmat(beforeAv_v,1,endThSprIdx-startThSprIdx+1);
            end
            
            if ~isempty(throatSpray_temp)
                throatSpraySuppressed = true;
            else
                throatSpraySuppressed = false;
            end
        end
        
        
        datatimes_proc = datatimes(procFirstIdx:procEndIdx-1);
        procedureDuration = minutes(datatimes(procEndIdx_raw) - datatimes(procStartIdx_raw));
        procedure = dataSource(:,procFirstIdx:procEndIdx-1);
        procedure_v = dataSource_v(:, procFirstIdx:procEndIdx-1);
        [procedure_max, maxIdx] = max(procedure,[],2);

        [procedure_max_sum, maxIdx_sum] = max(nansum(procedure(validDiams,:),1),[],2);
        [procedure_max_sum_v, maxIdx_sum_v] = max(nansum(procedure_v(validDiams,:),1),[],2);
        
        procedureTot = sum(procedure(validDiams,:),2);
        procedureTot_sum = sum(nansum(procedure(validDiams,:),1),2);
        procedureTot_sum_v = sum(nansum(procedure_v(validDiams,:),1),2);
        
        procedureMean = mean(procedure,2)/avSampleTime;
        procedureMean_sum = mean(nansum(procedure(validDiams,:),1),2)/avSampleTime;
        procedureMean_sum_v = mean(nansum(procedure_v(validDiams,:),1),2)/avSampleTime;
        
        refTime = datatimes_proc(1);
        
        eventThresh = 90; %In seconds
       
        for k = 1:size(maxIdx,1)
            eventTimes_val = eventTimes < datatimes_proc(maxIdx(k));
            
            if (nnz(eventTimes_val) == 0)
                nearestEventNames{k} = 'not recorded';
                tDiff(k) = NaN;
            else

                nearestEventIdxes(k) = knnsearch(seconds(eventTimes(eventTimes_val) - refTime),seconds(datatimes_proc(maxIdx(k)) - refTime));

                temp = table2cell(eventNames(nearestEventIdxes(k),1));
                nearestEventNames{k} = temp{1};
                tDiff(k) = datatimes_proc(maxIdx(k)) - eventTimes(nearestEventIdxes(k));
                
                if (tDiff(k) > seconds(eventThresh))
                    nearestEventNames{k} = 'not recorded';
                end
            end
        end
        
        eventTimes_val = eventTimes < datatimes_proc(maxIdx_sum);
        if (nnz(eventTimes_val) == 0)
            nearestEventNames_sum = 'not recorded';
            tDiff_sum = NaN;
        else
            nearestEventIdxes_sum = knnsearch(seconds(eventTimes(eventTimes_val) - refTime),seconds(datatimes_proc(maxIdx_sum) - refTime));
            temp = table2cell(eventNames(nearestEventIdxes_sum(1),1));
            nearestEventNames_sum = temp{1};
            tDiff_sum = datatimes_proc(maxIdx_sum) - eventTimes(nearestEventIdxes_sum);
            if (tDiff_sum > seconds(eventThresh))
                nearestEventNames_sum = 'not recorded';
            end
            
            
            if ~isempty(regexpi(nearestEventNames_sum, 'Null*.'))
                a = 1;
            end
        end
           
        eventTimes_val = eventTimes < datatimes_proc(maxIdx_sum_v);
        if (nnz(eventTimes_val) == 0)
            nearestEventNames_sum_v = 'not recorded';
            tDiff_sum_v = NaN;
        else
            nearestEventIdxes_sum_v = knnsearch(seconds(eventTimes(eventTimes_val) - refTime),seconds(datatimes_proc(maxIdx_sum_v) - refTime));
            temp = table2cell(eventNames(nearestEventIdxes_sum_v(1),1));
            nearestEventNames_sum_v = temp{1};
            tDiff_sum_v = datatimes_proc(maxIdx_sum_v) - eventTimes(nearestEventIdxes_sum_v);
            if (tDiff_sum_v > seconds(eventThresh))
                nearestEventNames_sum_v = 'not recorded';
            end
        end
        
        tempResults = table;
        tempResults.preProcedureMax_all = {preProcedure_max/avSampleTime};
        tempResults.preProcedureMax_sum = preProcedure_max_sum/avSampleTime;
        tempResults.preProcedureMax_sum_v = preProcedure_max_sum_v/avSampleTime;
        tempResults.preProcedureAv_sum = preProcedure_av/avSampleTime;
        
        tempResults.procedureMax_all = {procedure_max/avSampleTime};
        tempResults.procedureMax_sum = procedure_max_sum/avSampleTime;
        tempResults.procedureMax_sum_v = procedure_max_sum_v/avSampleTime;
        
        tempResults.procedureTot_all = {procedureTot};
        tempResults.procedureTot_sum = procedureTot_sum;
        tempResults.procedureTot_sum_v = procedureTot_sum_v;
        
        tempResults.procedureMean_all = {procedureMean};
        tempResults.procedureMean_sum = procedureMean_sum;
        tempResults.procedureMean_sum_v = procedureMean_sum_v;
        
        tempResults.nearestEvent_all = {nearestEventNames};
        tempResults.nearestEvent_sum = {nearestEventNames_sum};
        tempResults.nearestEvent_sum_v = {nearestEventNames_sum_v};
        
        tempResults.tDiff_all = {tDiff};
        tempResults.tDiff_sum = tDiff_sum;
        tempResults.tDiff_sum_v = tDiff_sum;
        
        tempResults.diameters = {diameters};
        tempResults.procedureDuration = procedureDuration;
        
        tempResults.StudyNumber = otherVars.StudyNumber;
        
        if suppressPosChanges
            tempResults.posChangesSuppressed = posChangesSuppressed;
        end
        
        if suppressThroatSpray
            tempResults.throatSpraySuppressed = throatSpraySuppressed;
        end
        
        tComb = join(otherVars,tempResults);
        
        tComb.procedureGroup = procedureGroup;
        if ~isempty(data_next)
            consecProcedure = true;
        else
            procedureGroup = procedureGroup+1;
        end
        
        % Test the BMI-years hypothesis
        tComb.BMIyears = tComb.BMI * tComb.Age;
        
        if isempty(procedureTable)
            procedureTable = tComb;
        else
            procedureTable = [procedureTable; tComb];
        end
        
        if ~isempty(data_next)
            tempResults_next = table;
            tempResults_next.data = {data_next};
            tempResults_next.time = {opTime2_next};
            tempResults_next.diameters = {diameters};
            
            dTime = opTime2_next - opTime2_next(1);
            idx_5min = find(dTime > minutes(5));
            
            validInterval = true;
            if ~(isempty(idx_5min))
                idx_5min = idx_5min(1);
            else
                validInterval = false;
            end
            idx_10min = find(dTime > minutes(10));
            if ~(isempty(idx_10min))
                idx_10min = idx_10min(1);
            else
                validInterval = false;
            end
            idx_20min = find(dTime > minutes(20));
            if ~(isempty(idx_20min))
                idx_20min = idx_20min(1);
            else
                validInterval = false;
            end
            
            if (validInterval)
                
                [bg_next, fg_next] = splitBGFG(data_next, avSampleTime, true(size(data_next,2),1));
                bg_next(bg_next < 0) = 0; %remove zeros
                
                avLevel = 5;
                
                idx_5min_end = idx_5min+avLevel-1;
                count_5min_all = mean(bg_next(:,idx_5min:idx_5min_end),2);%/avSampleTime;
                count_5min_sum = mean(sum(bg_next(validDiams,idx_5min:idx_5min_end),1),2);%/avSampleTime,1);

                idx_10min_end = idx_10min+avLevel-1;
                count_10min_all = mean(bg_next(:,idx_10min:idx_10min_end),2);%/avSampleTime;
                count_10min_sum = mean(sum(bg_next(validDiams,idx_10min:idx_10min_end),1),2);%/avSampleTime,1);

                idx_20min_end = idx_20min+avLevel-1;
                idx_20min_end = min(idx_20min_end, size(bg_next,2));
                count_20min_all = mean(bg_next(:,idx_20min:idx_20min_end),2);%/avSampleTime;
                count_20min_sum = mean(sum(bg_next(validDiams,idx_20min:idx_20min_end),1),2);%/avSampleTime,1);

%                 tempResults_next.count_5min_all = {count_5min_all./(bg_next(:,1)/avSampleTime)};
%                 tempResults_next.count_5min_sum = count_5min_sum./(nansum(bg_next(validDiams,1))/avSampleTime);
%                 tempResults_next.count_10min_all = {count_10min_all./(bg_next(:,1)/avSampleTime)};
%                 tempResults_next.count_10min_sum = count_10min_sum./(nansum(bg_next(validDiams,1))/avSampleTime);
%                 tempResults_next.count_20min_all = {count_20min_all./(bg_next(:,1)/avSampleTime)};
%                 tempResults_next.count_20min_sum = count_20min_sum./(nansum(bg_next(validDiams,1))/avSampleTime);
                
                tempResults_next.count_5min_all = {count_5min_all};
                tempResults_next.count_5min_sum = count_5min_sum;
                tempResults_next.count_10min_all = {count_10min_all};
                tempResults_next.count_10min_sum = count_10min_sum;
                tempResults_next.count_20min_all = {count_20min_all};
                tempResults_next.count_20min_sum = count_20min_sum;
                tempResults_next.procedureDuration = procedureDuration;

                %TODO repeat for each 
                logvals_norm = log(nansum(bg_next(validDiams,:),1)./nansum(bg_next(validDiams,1)));
                logvals = log(nansum(bg_next(validDiams,:),1)); %Do we need to divide by average sample time?
                
                logdiff = diff(logvals)/avSampleTime;
                logdiff_norm = diff(logvals_norm)/avSampleTime;
                validslope = logdiff < 0;
                validslope = [validslope(1), validslope];
                tempFiltBySlope = logvals;
                tempFiltBySlope(~validslope) = NaN;
                
                tempFiltBySlope_norm = logvals_norm;
                tempFiltBySlope_norm(~validslope) = NaN;
                
                tempSlope = [];
                lengths = [];
                slopes = [];
                count = 0;
                sizethresh = 3; % Max length of segment to fit slope to
                for kk = 1:size(tempFiltBySlope,2)
                    if isnan(tempFiltBySlope(kk)) || isinf(tempFiltBySlope(kk)) || count >= sizethresh
                        if count > 1
                            x = (0:size(tempSlope,1)-1).';
                            x = x*avSampleTime;
                            p = polyfit(x, tempSlope,1);
                            
                            if isnan(p(1))
                                a=1;
                            end
                            slopes = [slopes; p(1)];
                            lengths = [lengths; size(tempSlope,1)];
                        end
                        tempSlope = [];
                        count = 0;
                    end
                    
                    if ~isnan(tempFiltBySlope(kk)) && ~isinf(tempFiltBySlope(kk)) 
                        tempSlope = [tempSlope; tempFiltBySlope(kk)];
                        count = count+1;
                    end
                end
                
                tempResults_next.slopes = {slopes};
                tempResults_next.lengths = {lengths};
                
                if exist('slopes_temp')
                    slopes_temp = [slopes_temp; slopes];
                else
                    slopes_temp = slopes;
                end

                tempResults_next.StudyNumber = otherVars.StudyNumber;
                tComb_next = join(otherVars,tempResults_next);

                if isempty(interProcedureTable)
                    interProcedureTable = tComb_next;
                else
                    interProcedureTable = [interProcedureTable; tComb_next];
                end
                
                plotInterProc = false;
                
                if plotInterProc
                    subplot(3,2,1);
                    plot(minutes(opTime2_next - opTime2_next(1)), logvals_norm);
                    xlim([0,30]);
                    xlabel('Minutes after procedure end');
                    title('Raw data');
                    ylabel('No particles (log scale)');
                    hold on;

                    subplot(3,2,3);

                    plot(minutes(opTime2_next(1:end-1) - opTime2_next(1)), logdiff_norm);
                    xlim([0,30]);
                    xlabel('Minutes after procedure end');
                    hold on;

                    subplot(3,2,5)
                    plot(minutes(opTime2_next - opTime2_next(1)), tempFiltBySlope_norm);
                    xlim([0,30]);
                    xlabel('Minutes after procedure end');
                    title('Negative slopes only');
                    ylabel('No particles (log scale)');
                    hold on;

                    subplot(3,2,2)
                    histogram(slopes_temp,linspace(-0.01,0,30));
                    xlabel('Slopes');
                    pause(0.1);
                end
                
                
            end
            
        end

        % To do: store average results
        % To do: plot inter-procedure curves, 5 mins in, 10 mins in etc.
        
    end
    
end

upperGITable = procedureTable(procedureTable.Procedure == 'upper GI', :);
upperGITable_oral = procedureTable(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'oral', :);
upperGITable_nasal = procedureTable(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'nasal', :);
lowerGITable = procedureTable(procedureTable.Procedure == 'lower GI', :);
lowerGITable_colonoscopy = procedureTable(procedureTable.Procedure == 'lower GI' & procedureTable.ProcedureType == 'colonoscopy', :);
lowerGITable_sigmoid = procedureTable(procedureTable.Procedure == 'lower GI' & procedureTable.ProcedureType == 'sigmoidoscopy', :);

procGroups = categories(categorical(procedureTable.procedureGroup));
procGroupCounts = countcats(categorical(procedureTable.procedureGroup));
procGroupCounts_counts = countcats(categorical(procGroupCounts));
procGroupCounts_cats = categories(categorical(procGroupCounts));
bar(categorical(procGroupCounts_cats), (procGroupCounts_counts.*str2num(cell2mat(procGroupCounts_cats)))/size(procedureTable,1)*100);
xlabel('Procedure group size');
ylabel('Percentage of procedures');

%% General stats for DEN paper
age_ugi = upperGITable_oral.Age;
age_ugi_conv = age_ugi(upperGITable_oral.RoomType == 'endoscopy');
age_ugi_lam = age_ugi(upperGITable_oral.RoomType == 'laminar flow');
age_ugi_asent = age_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

[h_age_lam, p_age_lam, ci_age_lam] = ttest2(age_ugi_conv, age_ugi_lam);
disp(['UGI oral laminar flow vs conv age: p = ', num2str(p_age_lam)]);
[h_age_asent, p_age_asent, ci_age_asent] = ttest2(age_ugi_conv, age_ugi_asent);
disp(['UGI oral airsentry vs conv age: p = ', num2str(p_age_asent)]);

sex_ugi = upperGITable_oral.Sex;
sex_ugi_conv = sex_ugi(upperGITable_oral.RoomType == 'endoscopy');
sex_ugi_lam = sex_ugi(upperGITable_oral.RoomType == 'laminar flow');
sex_ugi_asent = sex_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

x = sum(sex_ugi_conv == 'female');
n = numel(sex_ugi_conv);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sex_conv_distp(p_idx) = binopdf(x,n,current_p);
end
sex_conv_distp = sex_conv_distp ./ (sum(sex_conv_distp)*dp);

x = sum(sex_ugi_lam == 'female');
n = numel(sex_ugi_lam);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sex_lam_distp(p_idx) = binopdf(x,n,current_p);
end
sex_lam_distp = sex_lam_distp ./ (sum(sex_lam_distp).*dp);
olap = min([sex_conv_distp;sex_lam_distp;]);
p_sex_lam = sum(olap.*dp);

% plot(p_vals, sex_conv_distp);
% hold on;
% plot(p_vals, sex_lam_distp, 'r');
% plot(p_vals, olap,'g');

x = sum(sex_ugi_asent == 'female');
n = numel(sex_ugi_asent);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sex_asent_distp(p_idx) = binopdf(x,n,current_p);
end
sex_asent_distp = sex_asent_distp ./ (sum(sex_asent_distp).*dp);
olap = min([sex_conv_distp;sex_asent_distp;]);
p_sex_asent = sum(olap.*dp);

disp(['UGI oral laminar flow vs conv sex : p = ', num2str(p_sex_lam)]);
disp(['UGI oral airsentry vs conv sex: p = ', num2str(p_sex_asent)]);


bmi_ugi = upperGITable_oral.BMI;
bmi_ugi_conv = bmi_ugi(upperGITable_oral.RoomType == 'endoscopy');
bmi_ugi_lam = bmi_ugi(upperGITable_oral.RoomType == 'laminar flow');
bmi_ugi_asent = bmi_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

[h_bmi_lam, p_bmi_lam, ci_bmi_lam] = ttest2(bmi_ugi_conv, bmi_ugi_lam);
disp(['UGI oral laminar flow vs conv bmi: p = ', num2str(p_bmi_lam)]);
[h_bmi_asent, p_bmi_asent, ci_bmi_asent] = ttest2(bmi_ugi_conv, bmi_ugi_asent);
disp(['UGI oral airsentry vs conv bmi: p = ', num2str(p_bmi_asent)]);

smoking_ugi = upperGITable_oral.Smoker;
smoking_ugi_conv = smoking_ugi(upperGITable_oral.RoomType == 'endoscopy');
smoking_ugi_lam = smoking_ugi(upperGITable_oral.RoomType == 'laminar flow');
smoking_ugi_asent = smoking_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

x = sum(smoking_ugi_conv == 'yes');
n = numel(smoking_ugi_conv);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    smoking_conv_distp(p_idx) = binopdf(x,n,current_p);
end
smoking_conv_distp = smoking_conv_distp ./ (sum(smoking_conv_distp)*dp);

x = sum(smoking_ugi_lam == 'yes');
n = numel(smoking_ugi_lam);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    smoking_lam_distp(p_idx) = binopdf(x,n,current_p);
end
smoking_lam_distp = smoking_lam_distp ./ (sum(smoking_lam_distp).*dp);
olap = min([smoking_conv_distp;smoking_lam_distp;]);
p_smoking_lam = sum(olap.*dp);

% plot(p_vals, smoking_conv_distp);
% hold on;
% plot(p_vals, smoking_lam_distp, 'r');
% plot(p_vals, olap,'g');

x = sum(smoking_ugi_asent == 'yes');
n = numel(smoking_ugi_asent);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    smoking_asent_distp(p_idx) = binopdf(x,n,current_p);
end
smoking_asent_distp = smoking_asent_distp ./ (sum(smoking_asent_distp).*dp);
olap = min([smoking_conv_distp;smoking_asent_distp;]);
p_smoking_asent = sum(olap.*dp);

disp(['UGI oral laminar flow vs conv smoking : p = ', num2str(p_smoking_lam)]);
disp(['UGI oral airsentry vs conv smoking: p = ', num2str(p_smoking_asent)]);


hiatus_ugi = upperGITable_oral.HiatusHernia;
hiatus_ugi_conv = hiatus_ugi(upperGITable_oral.RoomType == 'endoscopy');
hiatus_ugi_lam = hiatus_ugi(upperGITable_oral.RoomType == 'laminar flow');
hiatus_ugi_asent = hiatus_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

x = sum(hiatus_ugi_conv == 'yes');
n = numel(hiatus_ugi_conv);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    hiatus_conv_distp(p_idx) = binopdf(x,n,current_p);
end
hiatus_conv_distp = hiatus_conv_distp ./ (sum(hiatus_conv_distp)*dp);

x = sum(hiatus_ugi_lam == 'yes');
n = numel(hiatus_ugi_lam);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    hiatus_lam_distp(p_idx) = binopdf(x,n,current_p);
end
hiatus_lam_distp = hiatus_lam_distp ./ (sum(hiatus_lam_distp).*dp);
olap = min([hiatus_conv_distp;hiatus_lam_distp;]);
p_hiatus_lam = sum(olap.*dp);

% plot(p_vals, hiatus_conv_distp);
% hold on;
% plot(p_vals, hiatus_lam_distp, 'r');
% plot(p_vals, olap,'g');

x = sum(hiatus_ugi_asent == 'yes');
n = numel(hiatus_ugi_asent);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    hiatus_asent_distp(p_idx) = binopdf(x,n,current_p);
end
hiatus_asent_distp = hiatus_asent_distp ./ (sum(hiatus_asent_distp).*dp);
olap = min([hiatus_conv_distp;hiatus_asent_distp;]);
p_hiatus_asent = sum(olap.*dp);

disp(['UGI oral laminar flow vs conv hiatus : p = ', num2str(p_hiatus_lam)]);
disp(['UGI oral airsentry vs conv hiatus: p = ', num2str(p_hiatus_asent)]);



sedation_ugi = upperGITable_oral.Sedation;
sedation_ugi_conv = sedation_ugi(upperGITable_oral.RoomType == 'endoscopy');
sedation_ugi_lam = sedation_ugi(upperGITable_oral.RoomType == 'laminar flow');
sedation_ugi_asent = sedation_ugi(upperGITable_oral.RoomType == 'endo+airsentry');

x = sum(sedation_ugi_conv == 'midazolam');
n = numel(sedation_ugi_conv);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sedation_conv_distp(p_idx) = binopdf(x,n,current_p);
end
sedation_conv_distp = sedation_conv_distp ./ (sum(sedation_conv_distp)*dp);

x = sum(sedation_ugi_lam == 'midazolam');
n = numel(sedation_ugi_lam);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sedation_lam_distp(p_idx) = binopdf(x,n,current_p);
end
sedation_lam_distp = sedation_lam_distp ./ (sum(sedation_lam_distp).*dp);
olap = min([sedation_conv_distp;sedation_lam_distp;]);
p_sedation_lam = sum(olap.*dp);

% plot(p_vals, sedation_conv_distp);
% hold on;
% plot(p_vals, sedation_lam_distp, 'r');
% plot(p_vals, olap,'g');

x = sum(sedation_ugi_asent == 'midazolam');
n = numel(sedation_ugi_asent);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sedation_asent_distp(p_idx) = binopdf(x,n,current_p);
end
sedation_asent_distp = sedation_asent_distp ./ (sum(sedation_asent_distp).*dp);
olap = min([sedation_conv_distp;sedation_asent_distp;]);
p_sedation_asent = sum(olap.*dp);

disp(['UGI oral laminar flow vs conv sedation : p = ', num2str(p_sedation_lam)]);
disp(['UGI oral airsentry vs conv sedation: p = ', num2str(p_sedation_asent)]);

%% General stats for cytosponge paper
age_ugi = upperGITable_oral.Age;
age_ugi_ogd = age_ugi(upperGITable_oral.ProcedureType == 'gastroscopy');
age_ugi_cyto = age_ugi(upperGITable_oral.ProcedureType == 'cytosponge');

[h_age_cyto, p_age_cyto, ci_age_cyto] = ttest2(age_ugi_ogd, age_ugi_cyto);
disp(['Age OGD: ', num2str(min(age_ugi_ogd)), ' - ', num2str(max(age_ugi_ogd)), ' median ', num2str(median(age_ugi_ogd))]);
disp(['Age Cyto: ', num2str(min(age_ugi_cyto)), ' - ', num2str(max(age_ugi_cyto)), ' median ', num2str(median(age_ugi_cyto))]);
disp(['UGI gastroscopy vs cytosponge age: p = ', num2str(p_age_cyto)]);

sex_ugi = upperGITable_oral.Sex;
sex_ugi_ogd = sex_ugi(upperGITable_oral.ProcedureType == 'gastroscopy');
sex_ugi_cyto = sex_ugi(upperGITable_oral.ProcedureType == 'cytosponge');

x = sum(sex_ugi_ogd == 'female');
n = numel(sex_ugi_ogd);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sex_ogd_distp(p_idx) = binopdf(x,n,current_p);
end
sex_ogd_distp = sex_ogd_distp ./ (sum(sex_ogd_distp)*dp);

x = sum(sex_ugi_cyto == 'female');
n = numel(sex_ugi_cyto);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sex_cyto_distp(p_idx) = binopdf(x,n,current_p);
end
sex_cyto_distp = sex_cyto_distp ./ (sum(sex_cyto_distp).*dp);
olap = min([sex_ogd_distp;sex_cyto_distp;]);
p_sex_cyto = sum(olap.*dp);

disp(['Sex OGD: M - ', num2str(sum(sex_ugi_ogd == 'male')), ', F - ', num2str(sum(sex_ugi_ogd == 'female'))]);
disp(['Sex cyto: M - ', num2str(sum(sex_ugi_cyto == 'male')), ', F - ', num2str(sum(sex_ugi_cyto == 'female'))]);
disp(['UGI oral cytosonge vs ogd sex : p = ', num2str(p_sex_cyto)]);


bmi_ugi = upperGITable_oral.BMI;
bmi_ugi_ogd = bmi_ugi(upperGITable_oral.ProcedureType == 'gastroscopy');
bmi_ugi_cyto = bmi_ugi(upperGITable_oral.ProcedureType == 'cytosponge');

[h_bmi_cyto, p_bmi_cyto, ci_bmi_cyto] = ttest2(bmi_ugi_ogd, bmi_ugi_cyto);
disp(['BMI OGD: ', num2str(min(bmi_ugi_ogd)), ' - ', num2str(max(bmi_ugi_ogd)), ' median ', num2str(median(bmi_ugi_ogd))]);
disp(['BMI Cyto: ', num2str(min(bmi_ugi_cyto)), ' - ', num2str(max(bmi_ugi_cyto)), ' median ', num2str(nanmedian(bmi_ugi_cyto))]);
disp(['UGI oral cytosponge vs ogd bmi: p = ', num2str(p_bmi_cyto)]);

smoking_ugi = upperGITable_oral.Smoker;
smoking_ugi_ogd = smoking_ugi(upperGITable_oral.ProcedureType == 'gastroscopy');
smoking_ugi_cyto = smoking_ugi(upperGITable_oral.ProcedureType == 'cytosponge');

x = sum(smoking_ugi_ogd == 'yes');
n = numel(smoking_ugi_ogd);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    smoking_ogd_distp(p_idx) = binopdf(x,n,current_p);
end
smoking_ogd_distp = smoking_ogd_distp ./ (sum(smoking_ogd_distp)*dp);

x = sum(smoking_ugi_cyto == 'yes');
n = numel(smoking_ugi_cyto);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    smoking_cyto_distp(p_idx) = binopdf(x,n,current_p);
end
smoking_cyto_distp = smoking_cyto_distp ./ (sum(smoking_cyto_distp).*dp);
olap = min([smoking_ogd_distp;smoking_cyto_distp;]);
p_smoking_cyto = sum(olap.*dp);
disp(['Smoking OGD: Y - ', num2str(sum(smoking_ugi_ogd == 'yes')), ', no - ', num2str(sum(smoking_ugi_ogd == 'no'))]);
disp(['Smoking cyto: Y - ', num2str(sum(smoking_ugi_cyto == 'yes')), ', no - ', num2str(sum(smoking_ugi_cyto == 'no'))]);
disp(['UGI oral cytosponge vs ogd smoking : p = ', num2str(p_smoking_cyto)]);



sedation_ugi = upperGITable_oral.Sedation;
sedation_ugi_ogd = sedation_ugi(upperGITable_oral.ProcedureType == 'gastroscopy');
sedation_ugi_cyto = sedation_ugi(upperGITable_oral.ProcedureType == 'cytosponge');

%x = sum(sedation_ugi_ogd == 'throatspray');
x = numel(sedation_ugi_ogd);
n = numel(sedation_ugi_ogd);
p_vals = linspace(0,1,100);
dp = p_vals(2) - p_vals(1);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sedation_ogd_distp(p_idx) = binopdf(x,n,current_p);
end
sedation_ogd_distp = sedation_ogd_distp ./ (sum(sedation_ogd_distp)*dp);

x = sum(sedation_ugi_cyto == 'throat spray');
n = numel(sedation_ugi_cyto);
for p_idx = 1:size(p_vals,2)
    current_p = p_vals(p_idx);
    sedation_cyto_distp(p_idx) = binopdf(x,n,current_p);
end
sedation_cyto_distp = sedation_cyto_distp ./ (sum(sedation_cyto_distp).*dp);
olap = min([sedation_ogd_distp;sedation_cyto_distp;]);
p_sedation_cyto = sum(olap.*dp);

disp(['TS cyto: Y - ', num2str(sum(sedation_ugi_cyto == 'throat spray')), ', no - ', num2str(sum(sedation_ugi_cyto == 'none'))]);
disp(['UGI oral cytosponge vs ogd sedation : p = ', num2str(p_sedation_cyto)]);

%% Coughs room type
coughCount = upperGITable_oral.CoughEvents;
coughCount = coughCount./ upperGITable_oral.procedureDuration*20;
coughCount_ugi_conv = coughCount(upperGITable_oral.RoomType == 'endoscopy');
coughCount_ugi_lam = coughCount(upperGITable_oral.RoomType == 'laminar flow');
coughCount_ugi_asent = coughCount(upperGITable_oral.RoomType == 'endo+airsentry');

disp(['Mean no coughs: conv: ', num2str(mean(coughCount_ugi_conv)), ' laminar: ', num2str(mean(coughCount_ugi_lam)), ' air sentry: ', num2str(mean(coughCount_ugi_asent))]);

if ~isempty(coughCount_ugi_conv) && ~isempty(coughCount_ugi_lam)
    pd_c = fitdist(coughCount_ugi_conv, 'Poisson');
    pd_g = fitdist(coughCount_ugi_lam, 'Poisson');

    objFun2 = @(x) objFun(x,pd_c, pd_g);
    pVal = fmincon(objFun2,0.05,[],[],[],[],[1e-6],[0.999]);

    disp(['Cough count p val (conv v. lam): ', num2str(pVal), ',  ratio = ', num2str(mean(pd_c)/mean(pd_g))]);
end

if ~isempty(coughCount_ugi_conv) && ~isempty(coughCount_ugi_asent)
    pd_c = fitdist(coughCount_ugi_conv, 'Poisson');
    pd_g = fitdist(coughCount_ugi_asent, 'Poisson');

    objFun2 = @(x) objFun(x,pd_c, pd_g);
    pVal = fmincon(objFun2,0.05,[],[],[],[],[1e-6],[0.999]);

    disp(['Cough count p val (conv v. asent): ', num2str(pVal), ',  ratio = ', num2str(mean(pd_c)/mean(pd_g))]);
end

%% Coughs cyto
coughCount = upperGITable_oral.CoughEvents;
coughCount = coughCount./ upperGITable_oral.procedureDuration*20;
coughCount_ugi_endo = coughCount(upperGITable_oral.ProcedureType == 'gastroscopy');
coughCount_ugi_cyto = coughCount(upperGITable_oral.ProcedureType == 'cytosponge');

disp(['Mean no coughs: endo: ', num2str(mean(coughCount_ugi_endo)), ' cyto: ', num2str(mean(coughCount_ugi_cyto))]);

if ~isempty(coughCount_ugi_endo) && ~isempty(coughCount_ugi_cyto)
    pd_e = fitdist(coughCount_ugi_endo, 'Poisson');
    pd_c = fitdist(coughCount_ugi_cyto, 'Poisson');

    objFun2 = @(x) objFun(x,pd_e, pd_c);
    pVal = fmincon(objFun2,0.05,[],[],[],[],[1e-6],[0.999]);

    disp(['Cough count p val (endo v. cyto): ', num2str(pVal), ',  ratio = ', num2str(mean(pd_e)/mean(pd_c))]);
end

%%

procedureDuration = upperGITable_oral.procedureDuration;
procedureDuration_ugi_conv = procedureDuration(upperGITable_oral.RoomType == 'endoscopy');
procedureDuration_ugi_lam = procedureDuration(upperGITable_oral.RoomType == 'laminar flow');
procedureDuration_ugi_asent = procedureDuration(upperGITable_oral.RoomType == 'endo+airsentry');

disp(['Means: conv: ',  num2str(median(procedureDuration_ugi_conv)), ' laminar: ', num2str(median(procedureDuration_ugi_lam)), ' airsentry: ',  num2str(median(procedureDuration_ugi_asent))]);
[h_dur_lam, p_dur_lam, ci_dur_lam] = ttest2(log(procedureDuration_ugi_conv), log(procedureDuration_ugi_lam));
disp(['UGI oral laminar flow vs conv dur: p = ', num2str(p_dur_lam)]);
[h_dur_asent, p_dur_asent, ci_dur_asent] = ttest2(log(procedureDuration_ugi_conv), log(procedureDuration_ugi_asent));
disp(['UGI oral airsentry vs conv dur: p = ', num2str(p_dur_asent)]);



procedureDuration = upperGITable_oral.procedureDuration;
procedureDuration_ugi_ogd = procedureDuration(upperGITable_oral.ProcedureType == 'gastroscopy');
procedureDuration_ugi_cyto = procedureDuration(upperGITable_oral.ProcedureType == 'cytosponge');

disp(['Means: ogd: ',  num2str(median(procedureDuration_ugi_ogd)), ' cyto: ', num2str(median(procedureDuration_ugi_cyto))]);
[h_dur_lam, p_dur_cyto, ci_dur_cyto] = ttest2(log(procedureDuration_ugi_ogd), log(procedureDuration_ugi_cyto));
disp(['UGI oral cyto vs ogd dur: p = ', num2str(p_dur_cyto)]);
%%
if (~isempty(upperGITable_oral))
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1);
    pie(categorical(upperGITable_oral.nearestEvent_sum));
    title('Upper GI oral max events');
end


if (~isempty(upperGITable_nasal))
    subplot(1,3,2);
    pie(categorical(upperGITable_nasal.nearestEvent_sum));
    title('Upper GI nasal max events');
end

if (~isempty(lowerGITable))
    subplot(1,3,3);
    pie(categorical(lowerGITable.nearestEvent_sum));
    title('Lower GI max events');
end



saveas(gcf,fullfile(resultsFolder,['fullproc_pie_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_pie_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1);

if (~isempty(upperGITable_oral))
    pie(categorical(upperGITable_oral.nearestEvent_sum_v));
    title('Upper GI oral max events (volume)');
end

if (~isempty(upperGITable_nasal))
    subplot(1,3,2);
    pie(categorical(upperGITable_nasal.nearestEvent_sum_v));
    title('Upper GI nasal max events (volume)');
end

if (~isempty(lowerGITable))
    subplot(1,3,3);
    pie(categorical(lowerGITable.nearestEvent_sum_v));
    title('Lower GI max events (volume)');
end


saveas(gcf,fullfile(resultsFolder,['fullproc_pie_v_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_pie_v_', label, '.png']));

%%%
totPart = procedureTable.procedureTot_sum;
totPart_cat = procedureTable.Procedure;
route_cat = procedureTable.UGIroute;
type_cat = procedureTable.ProcedureType;

totPart_pre = procedureTable.preProcedureAv_sum .* procedureTable.procedureDuration * 60;
totPart_pre_cat = categorical(ones(size(totPart_pre)),[0,1],{'a','pre-procedure'});

[~, p_lp, ci_lp] = ttest(log(totPart(totPart_cat == 'lower GI')), log(totPart_pre(totPart_cat == 'lower GI')));
[~, p_uop, ci_uop] = ttest(log(totPart(totPart_cat == 'upper GI' & route_cat == 'oral')), log(totPart_pre(totPart_cat == 'upper GI'& route_cat == 'oral')));
[~, p_unp, ci_unp] = ttest(log(totPart(totPart_cat == 'upper GI' & route_cat == 'nasal')), log(totPart_pre(totPart_cat == 'upper GI'& route_cat == 'nasal')));
[~, p_luo, ci_luo] = ttest2(log(totPart(totPart_cat == 'lower GI')./procedureTable.procedureDuration(totPart_cat == 'lower GI')), log(totPart(totPart_cat == 'upper GI' & route_cat == 'oral')./procedureTable.procedureDuration(totPart_cat == 'upper GI'& route_cat == 'oral')));
[~, p_lun, ci_lun] = ttest2(log(totPart(totPart_cat == 'lower GI')./procedureTable.procedureDuration(totPart_cat == 'lower GI')), log(totPart(totPart_cat == 'upper GI' & route_cat == 'nasal')./procedureTable.procedureDuration(totPart_cat == 'upper GI'& route_cat == 'nasal')));
[~, p_uno, ci_uno] = ttest2(log(totPart(totPart_cat == 'upper GI' & route_cat == 'oral')./procedureTable.procedureDuration(totPart_cat == 'upper GI' & route_cat == 'oral')), log(totPart(totPart_cat == 'upper GI' & route_cat == 'nasal')./procedureTable.procedureDuration(totPart_cat == 'upper GI'& route_cat == 'nasal')));

[~, p_lcls, ci_lcls] = ttest2(log(totPart(totPart_cat == 'lower GI' & type_cat == 'colonoscopy')./procedureTable.procedureDuration(totPart_cat == 'lower GI' & type_cat == 'colonoscopy')), log(totPart(totPart_cat == 'lower GI' & type_cat == 'sigmoidoscopy')./procedureTable.procedureDuration(totPart_cat == 'lower GI'& type_cat == 'sigmoidoscopy')));



ratio_lp = exp(ci_lp);
ratio_uop = exp(ci_uop);
ratio_unp = exp(ci_unp);
ratio_luo = exp(ci_luo);
ratio_lun = exp(ci_lun);
ratio_uno = exp(ci_uno);
ratio_lcls = exp(ci_lcls);

meanLGI = exp(mean(log(totPart(totPart_cat == 'lower GI'))));
meanUGI_o = exp(mean(log(totPart(totPart_cat == 'upper GI' & route_cat == 'oral'))));
meanUGI_n = exp(mean(log(totPart(totPart_cat == 'upper GI' & route_cat == 'nasal'))));
meanLGI_colon = exp(mean(log(totPart(totPart_cat == 'lower GI' & type_cat == 'colonoscopy'))));
meanLGI_sigmoid = exp(mean(log(totPart(totPart_cat == 'lower GI' & type_cat == 'sigmoidoscopy'))));

mean_lp = sqrt(ratio_lp(1) * ratio_lp(2));
mean_uop = sqrt(ratio_uop(1) * ratio_uop(2));
mean_unp = sqrt(ratio_unp(1) * ratio_unp(2));
mean_luo = sqrt(ratio_luo(1) * ratio_luo(2));
mean_lun = sqrt(ratio_lun(1) * ratio_lun(2));
mean_uno = sqrt(ratio_uno(1) * ratio_uno(2));
mean_lcls = sqrt(ratio_lcls(1) * ratio_lcls(2));

disp(['Total LGI-pre ratio: ', num2str(mean_lp), ' (', num2str(ratio_lp(1)), '-', num2str(ratio_lp(2)), ') p=', num2str(p_lp), ', mean = ', num2str(meanLGI)]);
disp(['Total UGI-oral-pre ratio: ', num2str(mean_uop), ' (', num2str(ratio_uop(1)), '-', num2str(ratio_uop(2)), ') p=', num2str(p_uop), ', mean = ', num2str(meanUGI_o)]);
disp(['Total UGI-nasal-pre ratio: ', num2str(mean_unp), ' (', num2str(ratio_unp(1)), '-', num2str(ratio_unp(2)), ') p=', num2str(p_unp), ', mean = ', num2str(meanUGI_n)]);
disp(['Total LGI-UGI_o ratio: ', num2str(mean_luo), ' (', num2str(ratio_luo(1)), '-', num2str(ratio_luo(2)), ') p=', num2str(p_luo)]);
disp(['Colonoscopy mean: ', num2str(meanLGI_colon)]);
disp(['Sigmoidoscopy mean: ', num2str(meanLGI_sigmoid)]);

maxPart = procedureTable.procedureMax_sum;
maxPart_cat = procedureTable.Procedure;

maxPartPre = procedureTable.preProcedureMax_sum;
maxPartPre_cat = categorical(ones(size(maxPartPre)),[0,1],{'a','pre-procedure'});

maxPartPre = maxPartPre(logical(validPreEvent));
maxPartPre_cat = maxPartPre_cat(logical(validPreEvent));

maxPart_diff = maxPart(logical(validPreEvent)) - maxPartPre;
maxPart_diffCat = renamecats(maxPart_cat, {'upper GI diff', 'lower GI diff', 'N/A'});

maxPart_diffCat = maxPart_diffCat(logical(validPreEvent));

maxPart_diff_LGI = maxPart_diff(maxPart_diffCat == 'lower GI diff');
maxPart_diff_UGI = maxPart_diff(maxPart_diffCat == 'upper GI diff');

disp(['UGI: ', num2str(nnz(maxPart_diff_UGI > 0)), '/', num2str(numel(maxPart_diff_UGI)), ' = ', num2str(nnz(maxPart_diff_UGI > 0)/numel(maxPart_diff_UGI)*100), '%']);
disp(['LGI: ', num2str(nnz(maxPart_diff_LGI > 0)), '/', num2str(numel(maxPart_diff_LGI)), ' = ', num2str(nnz(maxPart_diff_LGI > 0)/numel(maxPart_diff_LGI)*100), '%']);

% [~, p_lp, ci_lp] = ttest2(log(maxPartPre), log(maxPart(maxPart_cat == 'lower GI')));
% [~, p_up, ci_up] = ttest2(log(maxPartPre), log(maxPart(maxPart_cat == 'upper GI')));
% [~, p_lu, ci_lu] = ttest2(log(maxPart(maxPart_cat == 'lower GI')), log(maxPart(maxPart_cat == 'upper GI')));
% 
% ratio_lp = exp(ci_lp);
% ratio_up = exp(ci_up);
% ratio_lu = exp(ci_lu);
% 
% mean_lp = sqrt(ratio_lp(1) * ratio_lp(2));
% mean_up = sqrt(ratio_up(1) * ratio_up(2));
% mean_lu = sqrt(ratio_lu(1) * ratio_lu(2));

meanPart = procedureTable.procedureMean_sum;
meanPart_cat = procedureTable.Procedure;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1);
%boxplot(totPart, totPart_cat);
totPart_cat2 = addcats(totPart_cat, 'nasal');
totPart_cat2(route_cat == 'nasal') = categorical({'nasal'},categories(totPart_cat2));

%% To avoid infinite, a hack
totPart_pre(totPart_pre == 0) = min(totPart_pre(totPart_pre > 0));

violinplot([totPart./totPart_pre], [totPart_cat2]);
title('Total no. particles');

text(0.5, max(totPart./totPart_pre)*0.8,['Total LGI-pre ratio: ', num2str(mean_lp), ' (', num2str(ratio_lp(1)), '-', num2str(ratio_lp(2)), ') p=', num2str(p_lp), ', mean = ', num2str(meanLGI)]);
text(0.5, max(totPart./totPart_pre)*0.7,['Total UGI_oral-pre ratio: ', num2str(mean_uop), ' (', num2str(ratio_uop(1)), '-', num2str(ratio_uop(2)), ') p=', num2str(p_uop), ', mean = ', num2str(meanUGI_o)]);
text(0.5, max(totPart./totPart_pre)*0.6,['Total UGI_nasal-pre ratio: ', num2str(mean_unp), ' (', num2str(ratio_unp(1)), '-', num2str(ratio_unp(2)), ') p=', num2str(p_unp), ', mean = ', num2str(meanUGI_n)]);
text(0.5, max(totPart./totPart_pre)*0.5,['Total LGI-UGI_o ratio: ', num2str(mean_luo), ' (', num2str(ratio_luo(1)), '-', num2str(ratio_luo(2)), ') p=', num2str(p_luo)]);
text(0.5, max(totPart./totPart_pre)*0.4,['Total LGI-UGI_n ratio: ', num2str(mean_lun), ' (', num2str(ratio_lun(1)), '-', num2str(ratio_lun(2)), ') p=', num2str(p_lun)]);
text(0.5, max(totPart./totPart_pre)*0.3,['Total UGI_o-UGI_n ratio: ', num2str(mean_uno), ' (', num2str(ratio_uno(1)), '-', num2str(ratio_uno(2)), ') p=', num2str(p_uno)]);

disp(['Total LGI-pre ratio: ', num2str(mean_lp), ' (', num2str(ratio_lp(1)), '-', num2str(ratio_lp(2)), ') p=', num2str(p_lp), ', mean = ', num2str(meanLGI)]);
disp(['Total UGI_oral-pre ratio: ', num2str(mean_uop), ' (', num2str(ratio_uop(1)), '-', num2str(ratio_uop(2)), ') p=', num2str(p_uop), ', mean = ', num2str(meanUGI_o)]);
disp(['Total UGI_nasal-pre ratio: ', num2str(mean_unp), ' (', num2str(ratio_unp(1)), '-', num2str(ratio_unp(2)), ') p=', num2str(p_unp), ', mean = ', num2str(meanUGI_n)]);
disp(['Total LGI-UGI_o ratio: ', num2str(mean_luo), ' (', num2str(ratio_luo(1)), '-', num2str(ratio_luo(2)), ') p=', num2str(p_luo)]);
disp(['Total LGI-UGI_n ratio: ', num2str(mean_lun), ' (', num2str(ratio_lun(1)), '-', num2str(ratio_lun(2)), ') p=', num2str(p_lun)]);
disp(['Total UGI_o-UGI_n ratio: ', num2str(mean_uno), ' (', num2str(ratio_uno(1)), '-', num2str(ratio_uno(2)), ') p=', num2str(p_uno)]);
disp(['Total LGI_c-LGI_s ratio: ', num2str(mean_lcls), ' (', num2str(ratio_lcls(1)), '-', num2str(ratio_lcls(2)), ') p=', num2str(p_lcls)]);

subplot(1,3,2);
%boxplot([maxPartPre; maxPart; maxPart_diff], [maxPartPre_cat; maxPart_cat; maxPart_diffCat]);
violinplot([maxPartPre; maxPart], [maxPartPre_cat; maxPart_cat]);
title('Max no. particles during proc./m^3s^{-1}');

subplot(1,3,3);
%boxplot(meanPart, meanPart_cat);
violinplot(meanPart, meanPart_cat);
title('Mean no. particles/m^3s^{-1}');

saveas(gcf,fullfile(resultsFolder,['fullproc_total_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_total_', label, '.png']));


%% Durations
durations = procedureTable.procedureDuration;
durations_LGI = durations(procedureTable.Procedure == 'lower GI');
durations_UGI = durations(procedureTable.Procedure == 'upper GI');

durations_colon = durations(procedureTable.Procedure == 'lower GI' & procedureTable.ProcedureType == 'colonoscopy');
durations_sigmoid = durations(procedureTable.Procedure == 'lower GI' & procedureTable.ProcedureType == 'sigmoidoscopy');

figure('units','normalized','outerposition',[0 0 1 1]);
histogram(durations_LGI,20);
hold on;
histogram(durations_UGI,20);
legend({'lower GI', 'upper GI'});
xlabel('procedure duration (minutes)');
ylabel('no. procedures');

mean_duration_LGI = mean(durations_LGI);
med_duration_LGI = median(durations_LGI);
mean_duration_UGI = mean(durations_UGI);
med_duration_UGI = median(durations_UGI);

mean_duration_colon = mean(durations_colon);
med_duration_colon = median(durations_colon);
mean_duration_sigmoid = mean(durations_sigmoid);
med_duration_sigmoid = median(durations_sigmoid);

disp(['UGI proc. duration: Mean: ', num2str(mean_duration_UGI), ', Med: ', num2str(med_duration_UGI), ' mins']);
disp(['LGI proc. duration: Mean: ', num2str(mean_duration_LGI), ', Med: ', num2str(med_duration_LGI), ' mins']);
disp(['Colonoscopy proc. duration: Mean: ', num2str(mean_duration_colon), ', Med: ', num2str(med_duration_colon), ' mins']);
disp(['Sigmoid proc. duration: Mean: ', num2str(mean_duration_sigmoid), ', Med: ', num2str(med_duration_sigmoid), ' mins']);

saveas(gcf,fullfile(resultsFolder,['fullproc_duration_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_duration_', label, '.png']));


figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
scatter(lowerGITable.procedureDuration,lowerGITable.procedureTot_sum);
temp_x = lowerGITable.procedureDuration;
temp_y = lowerGITable.procedureTot_sum;
valid = temp_y < 4e9;
temp_x = temp_x(valid);
temp_y = temp_y(valid);

if (~isempty(temp_x))
    [fitresult, gof, ~] = fit(temp_x,temp_y,'poly1');
    newx = linspace(min(lowerGITable.procedureDuration), max(lowerGITable.procedureDuration),100);
    yfit = feval(fitresult,newx);

    confidenceInt  = confint(fitresult);
    slopeconf = confidenceInt(:,1);

    if (size(lowerGITable.procedureTot_sum,1) > 2)
        p21 = predint(fitresult,newx,0.95,'functional','off');
    end
    hold on;
    plot(newx, yfit, 'k');
    
    
    if (size(lowerGITable.procedureTot_sum,1) > 2)
        plot(newx, p21, 'm--');
    end
    text(min(newx)*1.1, max(lowerGITable.procedureTot_sum)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);

    hold off;
    xlabel('procecure length (mins)');
    ylabel('tot no. particles');
    title('lower GI tot. particles v length');
end

subplot(1,2,2);
scatter(upperGITable_oral.procedureDuration,upperGITable_oral.procedureTot_sum);
[fitresult, gof, ~] = fit(upperGITable_oral.procedureDuration,upperGITable_oral.procedureTot_sum,'poly1');
newx = linspace(min(upperGITable_oral.procedureDuration), max(upperGITable_oral.procedureDuration),100);
yfit = feval(fitresult,newx);

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);

if (size(upperGITable_oral.procedureTot_sum,1) > 2)
    p21 = predint(fitresult,newx,0.95,'functional','off');
end
hold on;
plot(newx, yfit, 'k');

if (size(upperGITable_oral.procedureTot_sum,1) > 2)
    plot(newx, p21, 'm--');
end
text(min(newx)*1.1, max(upperGITable_oral.procedureTot_sum)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
hold off;
xlabel('procecure length (mins)');
ylabel('tot no. particles');
title('upper GI oral tot. particles v length');


saveas(gcf,fullfile(resultsFolder,['fullproc_durTrend_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_durTrend_', label, '.png']));


%% Plot vars

%% Anal tone
figure('units','normalized','outerposition',[0 0 1 1]);
analTone = procedureTable.AnalTone;
maxPart = procedureTable.procedureMax_sum;

analTone = analTone(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, analTone);
title('Max particles vs. Anal Tone');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(analTone == 'low')), log(maxPart(analTone == 'medium')));
[~, p_lh, ci_lh] = ttest2(log(maxPart(analTone == 'low')), log(maxPart(analTone == 'high')));
[~, p_mh, ci_mh] = ttest2(log(maxPart(analTone == 'medium')), log(maxPart(analTone == 'high')));

ci_lm = exp(ci_lm);
ci_lh = exp(ci_lh);
ci_mh = exp(ci_mh);

if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Low-Medium: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
    text(1.0, max(maxPart)*0.6,['Low-High: ', num2str(sqrt(ci_lh(1)*ci_lh(2))), ' (', num2str(ci_lh(1)), '-', num2str(ci_lh(2)), ') p = ', num2str(p_lh)]);
    text(1.0, max(maxPart)*0.4,['Medium-High: ', num2str(sqrt(ci_mh(1)*ci_mh(2))), ' (', num2str(ci_mh(1)), '-', num2str(ci_mh(2)), ') p = ', num2str(p_mh)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_analtone_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_analtone_', label, '.png']));


%% Sedation, proc. duration
figure('units','normalized','outerposition',[0 0 1 1]);
sedation = procedureTable.Sedation;
sedDuration = procedureTable.procedureDuration;

sedation = sedation(procedureTable.Procedure == 'upper GI');
sedDuration = sedDuration(procedureTable.Procedure == 'upper GI');
violinplot(sedDuration, sedation);
title('Duration particles vs. sedation, UGI');
ylabel('duration (mins)');

[~, p_sed, ci_sed] = ttest2(log(sedDuration(sedation == 'midazolam')), log(sedDuration(sedation == 'none')));
ci_sed = exp(ci_sed);

text(1.0, max(sedDuration)*0.8,['midazolam-none: ', num2str(sqrt(ci_sed(1)*ci_sed(2))), ' (', num2str(ci_sed(1)), '-', num2str(ci_sed(2)), ') p = ', num2str(p_sed)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_seddurugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_seddurugi_', label, '.png']));


figure('units','normalized','outerposition',[0 0 1 1]);
sedation = procedureTable.Sedation;
sedDuration = procedureTable.procedureDuration;

sedation = sedation(procedureTable.Procedure == 'lower GI');
sedDuration = sedDuration(procedureTable.Procedure == 'lower GI');
violinplot(sedDuration, sedation);
title('Duration particles vs. sedation, LGI');
ylabel('duration (mins)');

[~, p_sed, ci_sed] = ttest2(log(sedDuration(sedation == 'midazolam')), log(sedDuration(sedation == 'entonox')));
ci_sed = exp(ci_sed);

if (~isempty(sedDuration))
    text(1.0, max(sedDuration)*0.8,['midazolam-entonox: ', num2str(sqrt(ci_sed(1)*ci_sed(2))), ' (', num2str(ci_sed(1)), '-', num2str(ci_sed(2)), ') p = ', num2str(p_sed)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_seddurlgi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_seddurlgi_', label, '.png']));


%% Sex
figure('units','normalized','outerposition',[0 0 1 1]);
sex = procedureTable.Sex;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sex = sex(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, sex);
title('Max particles vs. sex upper GI');
ylabel('max no particles /m^3/s');

[~, p_mf, ci_mf] = ttest2(log(maxPart(sex == 'male')), log(maxPart(sex == 'female')));
ci_mf = exp(ci_mf);

text(1.0, max(maxPart)*0.8,['Male-female: ', num2str(sqrt(ci_mf(1)*ci_mf(2))), ' (', num2str(ci_mf(1)), '-', num2str(ci_mf(2)), ') p = ', num2str(p_mf)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_sexugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_sexugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
sex = procedureTable.Sex;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sex = sex(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, sex);
title('Max particles vs. sex lower GI');
ylabel('max no particles /m^3/s');

[~, p_mf, ci_mf] = ttest2(log(maxPart(sex == 'male')), log(maxPart(sex == 'female')));
ci_mf = 1./exp(ci_mf);

if ~isempty(maxPart)
    text(1.0, max(maxPart)*0.8,['Male-female: ', num2str(sqrt(ci_mf(1)*ci_mf(2))), ' (', num2str(ci_mf(1)), '-', num2str(ci_mf(2)), ') p = ', num2str(p_mf)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_sexlgi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_sexlgi_', label, '.png']));

%% Sedation
figure('units','normalized','outerposition',[0 0 1 1]);
sedation = procedureTable.Sedation;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sedation = sedation(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, sedation);
title('Max particles vs. sed upper');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(sedation == 'midazolam')), log(maxPart(sedation == 'none')));
ci_mn = exp(ci_mn);

[~, p_tsn, ci_tsn] = ttest2(log(maxPart(sedation == 'throat spray')), log(maxPart(sedation == 'none')));
ci_tsn = exp(ci_tsn);

[~, p_mts, ci_mts] = ttest2(log(maxPart(sedation == 'midazolam')), log(maxPart(sedation == 'throat spray')));
ci_mts = exp(ci_mts);

text(1.0, max(maxPart)*0.8,['Midaz-none: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);
text(1.0, max(maxPart)*0.6,['TS-none: ', num2str(sqrt(ci_tsn(1)*ci_tsn(2))), ' (', num2str(ci_tsn(1)), '-', num2str(ci_tsn(2)), ') p = ', num2str(p_tsn)]);
text(1.0, max(maxPart)*0.4,['Midaz-TS: ', num2str(sqrt(ci_mts(1)*ci_mts(2))), ' (', num2str(ci_mts(1)), '-', num2str(ci_mts(2)), ') p = ', num2str(p_mts)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_sedugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_sedugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
sedation = procedureTable.Sedation;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sedation = sedation(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, sedation);
title('Max particles vs. sed lower');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(sedation == 'midazolam')), log(maxPart(sedation == 'entonox')));

ci_mn = exp(ci_mn);

if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Midaz-entonox: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_sedlgi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_sedlgi_', label, '.png']));

%% Mask
figure('units','normalized','outerposition',[0 0 1 1]);
masked = procedureTable.PatientMask;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

masked = masked(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'oral');
maxPart = maxPart(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'oral');
violinplot(maxPart, masked);
title('Max particles vs. mask');
ylabel('max no particles /m^3/s');

[~, p_sb, ci_sb] = ttest2(log(maxPart(masked == 'surgical')), log(maxPart(masked == 'bronchoscopy')));
ci_sb = exp(ci_sb);
[~, p_sn, ci_sn] = ttest2(log(maxPart(masked == 'surgical')), log(maxPart(masked == 'no')));
ci_sn = exp(ci_sn);
[~, p_bn, ci_bn] = ttest2(log(maxPart(masked == 'bronchoscopy')), log(maxPart(masked == 'no')));
ci_bn = exp(ci_bn);

text(1.0, max(maxPart)*0.8,['Surgical-bronchoscopy: ', num2str(sqrt(ci_sb(1)*ci_sb(2))), ' (', num2str(ci_sb(1)), '-', num2str(ci_sb(2)), ') p = ', num2str(p_sb)]);
text(1.0, max(maxPart)*0.6,['Surgical-nonmasked: ', num2str(sqrt(ci_sn(1)*ci_sn(2))), ' (', num2str(ci_sn(1)), '-', num2str(ci_sn(2)), ') p = ', num2str(p_sn)]);
text(1.0, max(maxPart)*0.4,['Bronchoscopy-nonmasked: ', num2str(sqrt(ci_bn(1)*ci_bn(2))), ' (', num2str(ci_bn(1)), '-', num2str(ci_bn(2)), ') p = ', num2str(p_bn)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_maskedugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_maskedugi_', label, '.png']));

%% RoomType
figure('units','normalized','outerposition',[0 0 1 1]);
roomType = procedureTable.RoomType;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

roomType = roomType(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'oral');
maxPart = maxPart(procedureTable.Procedure == 'upper GI' & procedureTable.UGIroute == 'oral');
violinplot(maxPart, roomType);
title('Max particles vs. room type');
ylabel('max no particles /m^3/s');

endoData = maxPart(roomType == 'endoscopy');
lfData = maxPart(roomType == 'laminar flow');
asData = maxPart(roomType == 'endo+airsentry');

endo_mean = mean(log(endoData));
endo_std = std(log(endoData));

lfData_mean = mean(log(lfData));
asData_mean = mean(log(asData));

n_lf = sampsizepwr('t',[endo_mean, endo_std],lfData_mean);
pwrout_lft = sampsizepwr('t',[endo_mean, endo_std],lfData_mean,[],size(lfData,1));
disp(['To test laminar flow, with power 0.9, p=0.05 needs ', num2str(n_lf), '. Power with n=', num2str(size(lfData,1)), ' is ', num2str(pwrout_lft)]);


n_as = sampsizepwr('t',[endo_mean, endo_std],asData_mean);
pwrout_as = sampsizepwr('t',[endo_mean, endo_std],asData_mean,[],size(asData,1));
disp(['To test AS, with power 0.9, p=0.05 needs ', num2str(n_as), '. Power with n=', num2str(size(asData,1)), ' is ', num2str(pwrout_as)]);

[~, p_el, ci_el] = ttest2(log(maxPart(roomType == 'endoscopy')), log(maxPart(roomType == 'laminar flow')));
ci_el = exp(ci_el);
[~, p_ea, ci_ea] = ttest2(log(maxPart(roomType == 'endoscopy')), log(maxPart(roomType == 'endo+airsentry')));
ci_ea = exp(ci_ea);
[~, p_eal, ci_eal] = ttest2(log(maxPart(roomType == 'endo+airsentry')), log(maxPart(roomType == 'laminar flow')));
ci_eal = exp(ci_eal);

text(1.0, max(maxPart)*0.8,['Endo-laminar: ', num2str(sqrt(ci_el(1)*ci_el(2))), ' (', num2str(ci_el(1)), '-', num2str(ci_el(2)), ') p = ', num2str(p_el)]);
text(1.0, max(maxPart)*0.6,['Endo-endo-sent: ', num2str(sqrt(ci_ea(1)*ci_ea(2))), ' (', num2str(ci_ea(1)), '-', num2str(ci_ea(2)), ') p = ', num2str(p_ea)]);
text(1.0, max(maxPart)*0.4,['Endo-sent-laminar: ', num2str(sqrt(ci_eal(1)*ci_eal(2))), ' (', num2str(ci_eal(1)), '-', num2str(ci_eal(2)), ') p = ', num2str(p_eal)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_laminarugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_laminarugi_', label, '.png']));

%% Smoker
figure('units','normalized','outerposition',[0 0 1 1]);
smoker = procedureTable.Smoker;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

smoker = smoker(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, smoker);
title('Max particles vs. smoker');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(smoker == 'yes')), log(maxPart(smoker == 'no')));
ci_mn = exp(ci_mn);

text(1.0, max(maxPart)*0.8,['Smoker-nonsmoker: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_smokerugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_smokerugi_', label, '.png']));


figure('units','normalized','outerposition',[0 0 1 1]);
smoker = procedureTable.Smoker;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

smoker = smoker(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, smoker);
title('Max particles vs. smoker');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(smoker == 'yes')), log(maxPart(smoker == 'no')));

ci_mn = exp(ci_mn);
if ~isempty(maxPart)
    text(1.0, max(maxPart)*0.8,['Smoker-nonsmoker: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_smokerlgi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_smokerlgi_', label, '.png']));

%% Discomfort
figure('units','normalized','outerposition',[0 0 1 1]);
discomfort = procedureTable.Discomfort;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

discomfort = discomfort(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, discomfort);
title('Max particles vs. discom upper');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(discomfort == 'low')), log(maxPart(discomfort == 'medium')));
[~, p_lh, ci_lh] = ttest2(log(maxPart(discomfort == 'low')), log(maxPart(discomfort == 'high')));
[~, p_mh, ci_mh] = ttest2(log(maxPart(discomfort == 'medium')), log(maxPart(discomfort == 'high')));

ci_lm = exp(ci_lm);
ci_lh = exp(ci_lh);
ci_mh = exp(ci_mh);


text(1.0, max(maxPart)*0.8,['Low-Medium: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
text(1.0, max(maxPart)*0.6,['Low-High: ', num2str(sqrt(ci_lh(1)*ci_lh(2))), ' (', num2str(ci_lh(1)), '-', num2str(ci_lh(2)), ') p = ', num2str(p_lh)]);
text(1.0, max(maxPart)*0.4,['Medium-High: ', num2str(sqrt(ci_mh(1)*ci_mh(2))), ' (', num2str(ci_mh(1)), '-', num2str(ci_mh(2)), ') p = ', num2str(p_mh)]);

saveas(gcf,fullfile(resultsFolder,['fullproc_discomfortugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_discomfortugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
discomfort = procedureTable.Discomfort;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

discomfort = discomfort(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, discomfort);
title('Max particles vs. discom lower');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(discomfort == 'low')), log(maxPart(discomfort == 'medium')));
[~, p_lh, ci_lh] = ttest2(log(maxPart(discomfort == 'low')), log(maxPart(discomfort == 'high')));
[~, p_mh, ci_mh] = ttest2(log(maxPart(discomfort == 'medium')), log(maxPart(discomfort == 'high')));

ci_lm = exp(ci_lm);
ci_lh = exp(ci_lh);
ci_mh = exp(ci_mh);

if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Low-Medium: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
    text(1.0, max(maxPart)*0.6,['Low-High: ', num2str(sqrt(ci_lh(1)*ci_lh(2))), ' (', num2str(ci_lh(1)), '-', num2str(ci_lh(2)), ') p = ', num2str(p_lh)]);
    text(1.0, max(maxPart)*0.4,['Medium-High: ', num2str(sqrt(ci_mh(1)*ci_mh(2))), ' (', num2str(ci_mh(1)), '-', num2str(ci_mh(2)), ') p = ', num2str(p_mh)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_discomfortlgi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_discomfortlgi_', label, '.png']));
%% Co2
figure('units','normalized','outerposition',[0 0 1 1]);
CO2 = procedureTable.UseOfCO2orWater;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

CO2 = CO2(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, CO2);
title('Max particles vs. use of CO2');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(CO2 == 'CO2')), log(maxPart(CO2 == 'Water')));

ci_lm = exp(ci_lm);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['CO2-water: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_co2_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_co2_', label, '.png']));
%% Hysterectomy
figure('units','normalized','outerposition',[0 0 1 1]);
hyst = procedureTable.PreviousHysterectomy;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

hyst = hyst(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, hyst);
title('Max particles vs. previous hyst');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(hyst == 'yes')), log(maxPart(hyst == 'no')));

ci_lm = 1./exp(ci_lm);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Yes-No: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_hyst_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_hyst_', label, '.png']));
%% UGI route
figure('units','normalized','outerposition',[0 0 1 1]);
ugiroute = procedureTable.UGIroute;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

ugiroute = ugiroute(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, ugiroute);
title('Max particles vs. UGI route');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(ugiroute == 'nasal')), log(maxPart(ugiroute == 'oral')));
ci_lm = exp(ci_lm);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Nasal-Oral: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end

[~, p_lm, ci_lm] = ttest2(log(maxPart(ugiroute == 'nasal' | ugiroute == 'nasal abandoned')), log(maxPart(ugiroute == 'oral')));
ci_lm = exp(ci_lm);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.6,['Nasal (all) -Oral: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end
saveas(gcf,fullfile(resultsFolder,['fullproc_ugiroute_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_ugiroute_', label, '.png']));
%% UGI procedureType
figure('units','normalized','outerposition',[0 0 1 1]);
procedureType = procedureTable.ProcedureType;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

procedureType = procedureType(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, procedureType);
title('Max particles vs. procedureType');
ylabel('max no particles /m^3/s');

[~, p_gs, ci_gs] = ttest2(log(maxPart(procedureType == 'gastroscopy')), log(maxPart(procedureType == 'cytosponge')));
ci_gs = exp(ci_gs);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Gastroscopy-cytosponge: ', num2str(sqrt(ci_gs(1)*ci_gs(2))), ' (', num2str(ci_gs(1)), '-', num2str(ci_gs(2)), ') p = ', num2str(p_gs)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_proctype_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_proctype_', label, '.png']));

coughCount = upperGITable.CoughEvents;
coughCount = coughCount./ upperGITable.procedureDuration*20;
coughCount_c = coughCount(procedureType == 'cytosponge');

if ~cytospongeOnly
    if ~isempty(coughCount_c) && ~isempty(coughCount)
        pd_c = fitdist(coughCount_c, 'Poisson');
        coughCount_g = coughCount(procedureType == 'gastroscopy');
        pd_g = fitdist(coughCount_g, 'Poisson');

        objFun2 = @(x) objFun(x,pd_c, pd_g);
        pVal = fmincon(objFun2,0.05,[],[],[],[],[1e-6],[0.99999]);

        disp(['Cough count p val (cyto v. ogd): ', num2str(pVal), ',  ratio = ', num2str(mean(pd_c)/mean(pd_g))]);
    end
end

%% Hernia
figure('units','normalized','outerposition',[0 0 1 1]);
hiatus = upperGITable.HiatusHernia;
maxPart = upperGITable.procedureTot_sum ./ upperGITable.procedureDuration * 20;

violinplot(maxPart, hiatus);
title('Max particles vs. hiatus');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(hiatus == 'yes')), log(maxPart(hiatus == 'no')));

ci_lm = exp(ci_lm);
if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Yes-No: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_hernia_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_hernia_', label, '.png']));

coughCount = upperGITable.CoughEvents;
coughCount = coughCount./ upperGITable.procedureDuration*20;
coughCount_h = coughCount(hiatus == 'yes');

if ~isempty(coughCount_h) && ~isempty(coughCount)
    pd_h = fitdist(coughCount_h, 'Poisson');
    coughCount_nh = coughCount(hiatus == 'no');
    pd_nh = fitdist(coughCount_nh, 'Poisson');

    objFun2 = @(x) objFun(x,pd_h, pd_nh);
    pVal = fmincon(objFun2,0.05,[],[],[],[],[1e-6],[0.99999]);

    disp(['Cough count p val (hiatus v. non hiatus): ', num2str(pVal), ',  ratio = ', num2str(mean(pd_h)/mean(pd_nh))]);
end





%% Suctioning
figure('units','normalized','outerposition',[0 0 1 1]);
suctioning = procedureTable.Suctioning;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

suctioning = suctioning(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, suctioning);
title('Max particles vs. suctioning');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(suctioning == 'yes')), log(maxPart(suctioning == 'unknown')));

ci_lm = exp(ci_lm);

if (~isempty(maxPart))
    text(1.0, max(maxPart)*0.8,['Yes-no: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_suctioning_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_suctioning_', label, '.png']));
%% Age
figure('units','normalized','outerposition',[0 0 1 1]);
age = procedureTable.Age;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
age = age(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');

scatter(age, maxPart);
[fitresult, gof, ~] = fit(age(:),maxPart(:),'poly1');
newx = linspace(min(age), max(age),100);
yfit = feval(fitresult,newx);

if (size(maxPart,1) > 2)
    p21 = predint(fitresult,newx,0.95,'functional','off');
end
hold on;
plot(newx, yfit, 'k');

if (size(maxPart,1) > 2)
    plot(newx, p21, 'm--');
end
%text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
hold off;
xlabel('age');
ylabel('max no. particles');
title('age upper GI');

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);
if (~isempty(maxPart))
    text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_ageugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_ageugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
age = procedureTable.Age;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
age = age(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');

if (~isempty(maxPart))
    scatter(age, maxPart);
    [fitresult, gof, ~] = fit(age(:),maxPart(:),'poly1');
    newx = linspace(min(age), max(age),100);
    yfit = feval(fitresult,newx);

    if (size(maxPart,1) > 2)
        p21 = predint(fitresult,newx,0.95,'functional','off');
    end
    hold on;
    plot(newx, yfit, 'k');

    if (size(maxPart,1) > 2)
        plot(newx, p21, 'm--');
    end
    %text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
    hold off;
    xlabel('age');
    ylabel('max no. particles');
    title('age lower GI');

    confidenceInt  = confint(fitresult);
    slopeconf = confidenceInt(:,1);
    if (~isempty(maxPart))
        text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
    end
    saveas(gcf,fullfile(resultsFolder,['fullproc_agelgi_', label, '.fig']));
    saveas(gcf,fullfile(resultsFolder,['fullproc_agelgi_', label, '.png']));
end

%% BMI
figure('units','normalized','outerposition',[0 0 1 1]);
BMI = procedureTable.BMI;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
BMI = BMI(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');

val = ~isnan(BMI);
BMI = BMI(val);
maxPart = maxPart(val);

scatter(BMI, maxPart);
[fitresult, gof, ~] = fit(BMI(:),maxPart(:),'poly1');
newx = linspace(min(BMI), max(BMI),100);
yfit = feval(fitresult,newx);

if (size(maxPart,1) > 2)
    p21 = predint(fitresult,newx,0.95,'functional','off');
end
hold on;
plot(newx, yfit, 'k');

if (size(maxPart,1) > 2)
    plot(newx, p21, 'm--');
end
%text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
hold off;
xlabel('BMI');
ylabel('max no. particles');
title('BMI upper GI');

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);
if (~isempty(maxPart))
    text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_bmiugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_bmiugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
BMI = procedureTable.BMI;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
BMI = BMI(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');

if ~isempty(BMI)
    val = ~isnan(BMI);
    BMI = BMI(val);
    maxPart = maxPart(val);

    scatter(BMI, maxPart);
    [fitresult, gof, ~] = fit(BMI(:),maxPart(:),'poly1');
    newx = linspace(min(BMI), max(BMI),100);
    yfit = feval(fitresult,newx);

    if (size(maxPart,1) > 2)
        p21 = predint(fitresult,newx,0.95,'functional','off');
    end
    hold on;
    plot(newx, yfit, 'k');

    if (size(maxPart,1) > 2)
        plot(newx, p21, 'm--');
    end
    %text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
    hold off;
    xlabel('BMI');
    ylabel('max no. particles');
    title('BMI lower GI');
    a = 1;

    confidenceInt  = confint(fitresult);
    slopeconf = confidenceInt(:,1);
    if (~isempty(maxPart))
        text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
    end
    saveas(gcf,fullfile(resultsFolder,['fullproc_bmilgi_', label, '.fig']));
    saveas(gcf,fullfile(resultsFolder,['fullproc_bmilgi_', label, '.png']));
end

%% BMI years
figure('units','normalized','outerposition',[0 0 1 1]);
BMIyears = procedureTable.BMIyears;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
BMIyears = BMIyears(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');

val = ~isnan(BMIyears);
BMIyears = BMIyears(val);
maxPart = maxPart(val);

scatter(BMIyears, maxPart);
[fitresult, gof, ~] = fit(BMIyears(:),maxPart(:),'poly1');
newx = linspace(min(BMIyears), max(BMIyears),100);
yfit = feval(fitresult,newx);

if (size(maxPart,1) > 2)
    p21 = predint(fitresult,newx,0.95,'functional','off');
end
hold on;
plot(newx, yfit, 'k');

if (size(maxPart,1) > 2)
    plot(newx, p21, 'm--');
end
%text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
hold off;
xlabel('BMIyears');
ylabel('max no. particles');
title('BMIyears upper GI');

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);
if (~isempty(maxPart))
    text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
end

saveas(gcf,fullfile(resultsFolder,['fullproc_BMIyearsugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_BMIyearsugi_', label, '.png']));

figure('units','normalized','outerposition',[0 0 1 1]);
BMIyears = procedureTable.BMIyears;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
BMIyears = BMIyears(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');

if ~isempty(BMIyears)
    val = ~isnan(BMIyears);
    BMIyears = BMIyears(val);
    maxPart = maxPart(val);

    scatter(BMIyears, maxPart);
    [fitresult, gof, ~] = fit(BMIyears(:),maxPart(:),'poly1');
    newx = linspace(min(BMIyears), max(BMIyears),100);
    yfit = feval(fitresult,newx);

    if (size(maxPart,1) > 2)
        p21 = predint(fitresult,newx,0.95,'functional','off');
    end
    hold on;
    plot(newx, yfit, 'k');

    if (size(maxPart,1) > 2)
        plot(newx, p21, 'm--');
    end
    %text(min(newx)*1.1, max(maxPart)*0.8,['r = ', num2str(gof.rsquare)]);
    hold off;
    xlabel('BMIyears');
    ylabel('max no. particles');
    title('BMIyears lower GI');
    a = 1;

    confidenceInt  = confint(fitresult);
    slopeconf = confidenceInt(:,1);
    if (~isempty(maxPart))
        text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
    end

    saveas(gcf,fullfile(resultsFolder,['fullproc_BMIyearslgi_', label, '.fig']));
    saveas(gcf,fullfile(resultsFolder,['fullproc_BMIyearslgi_', label, '.png']));
end

%% Interprocedure
usesAirSentry = interProcedureTable.AirSentryUsed;

part5min = interProcedureTable(~usesAirSentry,:).count_5min_sum;
part5min_cat = interProcedureTable(~usesAirSentry,:).Procedure;

part10min = interProcedureTable(~usesAirSentry,:).count_10min_sum;
part10min_cat = interProcedureTable(~usesAirSentry,:).Procedure;

part20min = interProcedureTable(~usesAirSentry,:).count_20min_sum;
part20min_cat = interProcedureTable(~usesAirSentry,:).Procedure;

part5min_sent = interProcedureTable(usesAirSentry,:).count_5min_sum;
part5min_sent_cat = interProcedureTable(usesAirSentry,:).Procedure;

part10min_sent = interProcedureTable(usesAirSentry,:).count_10min_sum;
part10min_sent_cat = interProcedureTable(usesAirSentry,:).Procedure;

part20min_sent = interProcedureTable(usesAirSentry,:).count_20min_sum;
part20min_sent_cat = interProcedureTable(usesAirSentry,:).Procedure;

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
temp_data = [part5min; part10min; part20min];
temp_cats = [part5min_cat; part10min_cat; part20min_cat;];
temp_cats_new = [5*ones(size(part5min)), 10*ones(size(part10min)), 20*ones(size(part20min))];

temp_data_lg = temp_data(temp_cats == 'lower GI');
temp_data_lg_cats = temp_cats_new(temp_cats == 'lower GI');

violinplot(temp_data_lg, temp_data_lg_cats);
title('Lower GI');
ylimFirst = ylim;
ylim([0, max(ylimFirst)]);
xlabel('Mins after end of procedure');
ylabel('No. particles');

subplot(2,2,2);
temp_data = [part5min; part10min; part20min];
temp_cats = [part5min_cat; part10min_cat; part20min_cat;];
temp_cats_new = [5*ones(size(part5min)), 10*ones(size(part10min)), 20*ones(size(part20min))];

temp_data_ug = temp_data(temp_cats == 'upper GI');
temp_data_ug_cats = temp_cats_new(temp_cats == 'upper GI');

violinplot(temp_data_ug, temp_data_ug_cats);
title('Upper GI');
ylim(ylimFirst);
ylim([0, max(ylimFirst)]);
xlabel('Mins after end of procedure');
ylabel('No. particles');

subplot(2,2,3);
temp_data = [part5min_sent; part10min_sent; part20min_sent];
temp_cats = [part5min_sent_cat; part10min_sent_cat; part20min_sent_cat;];
temp_cats_new = [5*ones(size(part5min_sent)), 10*ones(size(part10min_sent)), 20*ones(size(part20min_sent))];

temp_data_lg = temp_data(temp_cats == 'lower GI');
temp_data_lg_cats = temp_cats_new(temp_cats == 'lower GI');

violinplot(temp_data_lg, temp_data_lg_cats);
title('Lower GI with AirSentry');
ylim(ylimFirst);
ylim([0, max(ylimFirst)]);
xlabel('Mins after end of procedure');
ylabel('No. particles');

subplot(2,2,4);
temp_data = [part5min_sent; part10min_sent; part20min_sent];
temp_cats = [part5min_sent_cat; part10min_sent_cat; part20min_sent_cat;];
temp_cats_new = [5*ones(size(part5min_sent)), 10*ones(size(part10min_sent)), 20*ones(size(part20min_sent))];

temp_data_ug = temp_data(temp_cats == 'upper GI');
temp_data_ug_cats = temp_cats_new(temp_cats == 'upper GI');

violinplot(temp_data_ug, temp_data_ug_cats);
title('Upper GI with AirSentry');
ylim(ylimFirst);
ylim([0, max(ylimFirst)]);
xlabel('Mins after end of procedure');
ylabel('No. particles');

saveas(gcf,fullfile(resultsFolder,['fullproc_interproctot_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_interproctot_', label, '.png']));

%% Now compute slopes (upper+lower)
medianDuration = median(interProcedureTable.procedureDuration);
iqrDuration = iqr(interProcedureTable.procedureDuration);
disp(['Preceeding procedure (upper+lower) has median = ', num2str(medianDuration), ' and IQR = ', num2str(iqrDuration)]);
slopes = interProcedureTable.slopes;
lengths = interProcedureTable.lengths;
slopes_nosent = slopes(~interProcedureTable.AirSentryUsed);
lengths_nosent = lengths(~interProcedureTable.AirSentryUsed);
slopes_nosent = vertcat(slopes_nosent{:});
lengths_nosent = vertcat(lengths_nosent{:});

slopes_nosent_f = [];
for k=1:size(slopes_nosent,1)
    slopes_nosent_f = [slopes_nosent_f; ones(lengths_nosent(k),1)*slopes_nosent(k)];
end

slopes_sent = slopes(interProcedureTable.AirSentryUsed);
lengths_sent = lengths(interProcedureTable.AirSentryUsed);
slopes_sent = vertcat(slopes_sent{:});
lengths_sent = vertcat(lengths_sent{:});

slopes_sent_f = [];
for k=1:size(slopes_sent,1)
    slopes_sent_f = [slopes_sent_f; ones(lengths_sent(k),1)*slopes_sent(k)];
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
histogram(slopes_nosent_f,50);
nosent_mean = sum(lengths_nosent/sum(lengths_nosent) .* slopes_nosent);
sent_mean = sum(lengths_sent/sum(lengths_sent) .* slopes_sent);
nosent_med = median(slopes_nosent);
sent_med = median(slopes_sent);
xlabel('slope');
ylabel('frequency');
title(['slopes w/o and w airsentry (upper and lower) (medians = ', num2str(nosent_med), ', ', num2str(sent_med), ')']);

hold on;
histogram(slopes_sent_f,20);

%xlabel('slope');
%ylabel('frequency');
%title(['slopes with airsentry (mean = ', num2str(sent_mean)]);
hold off;

subplot(1,3,2);
cats = [zeros(size(slopes_nosent,1),1);ones(size(slopes_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([slopes_nosent; slopes_sent], cats);
ylabel('Decay constant (s^{-1})');

% Compute confidence interval
if (~isempty(slopes_sent))
    pval_f = 0;
    ci_f = [];
    h_cum = 0;
    for k=1:1000
        temp_nosent = datasample(slopes_nosent,size(slopes_nosent,1),'Weights',lengths_nosent);
        temp_sent = datasample(slopes_sent,size(slopes_sent,1)*2,'Weights',lengths_sent);

        %[h,pval_temp, ci_temp] = ttest2(temp_nosent, temp_sent, 'tail', 'both');
        %[h,pval_temp, ci_temp] = kstest2(temp_nosent, temp_sent, 'Tail', 'unequal'); %Since they are highly non-normal (anderson-darling test) 
        [pval_temp, h, ci_temp] = ranksum(temp_nosent, temp_sent, 'tail', 'right'); %Since they are highly non-normal (anderson-darling test)
        ci_temp = [0,0];

        pval_f = pval_f + pval_temp;

        if isempty(ci_f)
            ci_f = ci_temp;
        else
            ci_f = ci_f + ci_temp;
        end
        
        h_cum = h_cum+h;
    end
    pval_f = pval_f/1000;
    ci_f = ci_f/1000;
    h_cum = h_cum/1000;
    
    disp(['Study power for slope method: ', num2str(h_cum)]);
end

partLevel = 0.5;
t_95_nosent = log(partLevel)./slopes_nosent_f/60;
t_95_sent = log(partLevel)./slopes_sent_f/60;

t_95_nosent_med = log(partLevel)/nosent_med/60;
t_95_sent_med = log(partLevel)/sent_med/60;
%t_95_sent_mean_lb = log(partLevel)/(nosent_mean - ci_f(1))/60;
%t_95_sent_mean_ub = log(partLevel)/(nosent_mean - ci_f(2))/60;

%ub2 = (nosent_mean - ci_f(1))/nosent_mean;
m2 = (sent_med)/nosent_med;
%lb2 = (nosent_mean - ci_f(2))/nosent_mean;



subplot(1,3,3);
cats = [zeros(size(t_95_nosent,1),1);ones(size(t_95_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([t_95_nosent; t_95_sent],cats);
ylabel('minutes');
title('time taken to reach 50%');
ylim([0,120]);
if (~isempty(slopes_sent))
    text(0,110,['p=',num2str(pval_f), ', med ratio=', num2str(m2), 't_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
    text(0,100,['n_{proc}=', num2str(nnz(usesAirSentry)), '/', num2str(nnz(~usesAirSentry)), ':  ', num2str(numel(lengths_sent)), '/', num2str(numel(lengths_nosent))]);
else
    text(0,110,['t_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
end
    
%histogram(t_95_nosent, 100);
%hold on;
%histogram(t_95_sent, 50);
%hold off;

saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_upper+lower_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_upper+lower', label, '.png']));

%%
%% Now compute slopes (lower)
interProcedureTable_lower = interProcedureTable(interProcedureTable.Procedure == 'lower GI',:);
usesAirSentry = interProcedureTable_lower.AirSentryUsed;
medianDuration = median(interProcedureTable_lower.procedureDuration);
iqrDuration = iqr(interProcedureTable_lower.procedureDuration);
disp(['Preceeding procedure (lower) has median = ', num2str(medianDuration), ' and IQR = ', num2str(iqrDuration)]);
slopes = interProcedureTable_lower.slopes;
lengths = interProcedureTable_lower.lengths;
slopes_nosent = slopes(~interProcedureTable_lower.AirSentryUsed);
lengths_nosent = lengths(~interProcedureTable_lower.AirSentryUsed);
slopes_nosent = vertcat(slopes_nosent{:});
lengths_nosent = vertcat(lengths_nosent{:});

slopes_nosent_f = [];
for k=1:size(slopes_nosent,1)
    slopes_nosent_f = [slopes_nosent_f; ones(lengths_nosent(k),1)*slopes_nosent(k)];
end

slopes_sent = slopes(interProcedureTable_lower.AirSentryUsed);
lengths_sent = lengths(interProcedureTable_lower.AirSentryUsed);
slopes_sent = vertcat(slopes_sent{:});
lengths_sent = vertcat(lengths_sent{:});

slopes_sent_f = [];
for k=1:size(slopes_sent,1)
    slopes_sent_f = [slopes_sent_f; ones(lengths_sent(k),1)*slopes_sent(k)];
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
histogram(slopes_nosent_f,50);
nosent_mean = sum(lengths_nosent/sum(lengths_nosent) .* slopes_nosent);
sent_mean = sum(lengths_sent/sum(lengths_sent) .* slopes_sent);
nosent_med = median(slopes_nosent);
sent_med = median(slopes_sent);
xlabel('slope');
ylabel('frequency');
title(['slopes w/o and w airsentry (lower) (medians = ', num2str(nosent_med), ', ', num2str(sent_med), ')']);

hold on;
histogram(slopes_sent_f,20);

%xlabel('slope');
%ylabel('frequency');
%title(['slopes with airsentry (mean = ', num2str(sent_mean)]);
hold off;

subplot(1,3,2);
cats = [zeros(size(slopes_nosent,1),1);ones(size(slopes_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([slopes_nosent; slopes_sent], cats);
ylabel('Decay constant (s^{-1})');

% Compute confidence interval
if (~isempty(slopes_sent))
    pval_f = 0;
    ci_f = [];
    for k=1:1000
        temp_nosent = datasample(slopes_nosent,size(slopes_nosent,1),'Weights',lengths_nosent);
        temp_sent = datasample(slopes_sent,size(slopes_sent,1)*2,'Weights',lengths_sent);

        %[h,pval_temp, ci_temp] = ttest2(temp_nosent, temp_sent, 'tail', 'both');
        %[h,pval_temp, ci_temp] = kstest2(temp_nosent, temp_sent, 'Tail', 'unequal'); %Since they are highly non-normal (anderson-darling test) 
        [pval_temp, h, ci_temp] = ranksum(temp_nosent, temp_sent, 'tail', 'right'); %Since they are highly non-normal (anderson-darling test)
        ci_temp = [0,0];

        pval_f = pval_f + pval_temp;

        if isempty(ci_f)
            ci_f = ci_temp;
        else
            ci_f = ci_f + ci_temp;
        end
    end
    pval_f = pval_f/1000;
    ci_f = ci_f/1000;



    [h,pval, ci] = ttest2(slopes_nosent, slopes_sent, 'tail', 'both');
end

partLevel = 0.5;
t_95_nosent = log(partLevel)./slopes_nosent_f/60;
t_95_sent = log(partLevel)./slopes_sent_f/60;

t_95_nosent_med = log(partLevel)/nosent_med/60;
t_95_sent_med = log(partLevel)/sent_med/60;
%t_95_sent_mean_lb = log(partLevel)/(nosent_mean - ci_f(1))/60;
%t_95_sent_mean_ub = log(partLevel)/(nosent_mean - ci_f(2))/60;

%ub2 = (nosent_mean - ci_f(1))/nosent_mean;
m2 = (sent_med)/nosent_med;
%lb2 = (nosent_mean - ci_f(2))/nosent_mean;



subplot(1,3,3);
cats = [zeros(size(t_95_nosent,1),1);ones(size(t_95_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([t_95_nosent; t_95_sent],cats);
ylabel('minutes');
title('time taken to reach 50%');
ylim([0,120]);
if (~isempty(slopes_sent))
    text(0,110,['p=',num2str(pval_f), ', med ratio=', num2str(m2), 't_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
    text(0,100,['n_{proc}=', num2str(nnz(usesAirSentry)), '/', num2str(nnz(~usesAirSentry)), ':  ', num2str(numel(lengths_sent)), '/', num2str(numel(lengths_nosent))]);
else
    text(0,110,['t_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
end
    
%histogram(t_95_nosent, 100);
%hold on;
%histogram(t_95_sent, 50);
%hold off;

saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_lower', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_lower', label, '.png']));
%%
%% Now compute slopes (upper)
interProcedureTable_upper = interProcedureTable(interProcedureTable.Procedure == 'upper GI',:);
usesAirSentry = interProcedureTable_upper.AirSentryUsed;
medianDuration = median(interProcedureTable_upper.procedureDuration);
iqrDuration = iqr(interProcedureTable_upper.procedureDuration);
disp(['Preceeding procedure (upper) has median = ', num2str(medianDuration), ' and IQR = ', num2str(iqrDuration)]);
slopes = interProcedureTable.slopes;
lengths = interProcedureTable.lengths;
slopes_nosent = slopes(~interProcedureTable.AirSentryUsed);
lengths_nosent = lengths(~interProcedureTable.AirSentryUsed);
slopes_nosent = vertcat(slopes_nosent{:});
lengths_nosent = vertcat(lengths_nosent{:});

slopes_nosent_f = [];
for k=1:size(slopes_nosent,1)
    slopes_nosent_f = [slopes_nosent_f; ones(lengths_nosent(k),1)*slopes_nosent(k)];
end

slopes_sent = slopes(interProcedureTable.AirSentryUsed);
lengths_sent = lengths(interProcedureTable.AirSentryUsed);
slopes_sent = vertcat(slopes_sent{:});
lengths_sent = vertcat(lengths_sent{:});

slopes_sent_f = [];
for k=1:size(slopes_sent,1)
    slopes_sent_f = [slopes_sent_f; ones(lengths_sent(k),1)*slopes_sent(k)];
end

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
histogram(slopes_nosent_f,50);
nosent_mean = sum(lengths_nosent/sum(lengths_nosent) .* slopes_nosent);
sent_mean = sum(lengths_sent/sum(lengths_sent) .* slopes_sent);
nosent_med = median(slopes_nosent);
sent_med = median(slopes_sent);
xlabel('slope');
ylabel('frequency');
title(['slopes w/o and w airsentry (upper) (medians = ', num2str(nosent_med), ', ', num2str(sent_med), ')']);

hold on;
histogram(slopes_sent_f,20);

%xlabel('slope');
%ylabel('frequency');
%title(['slopes with airsentry (mean = ', num2str(sent_mean)]);
hold off;

subplot(1,3,2);
cats = [zeros(size(slopes_nosent,1),1);ones(size(slopes_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([slopes_nosent; slopes_sent], cats);
ylabel('Decay constant (s^{-1})');

% Compute confidence interval
if (~isempty(slopes_sent))
    pval_f = 0;
    ci_f = [];
    for k=1:1000
        temp_nosent = datasample(slopes_nosent,size(slopes_nosent,1),'Weights',lengths_nosent);
        temp_sent = datasample(slopes_sent,size(slopes_sent,1)*2,'Weights',lengths_sent);

        %[h,pval_temp, ci_temp] = ttest2(temp_nosent, temp_sent, 'tail', 'both');
        %[h,pval_temp, ci_temp] = kstest2(temp_nosent, temp_sent, 'Tail', 'unequal'); %Since they are highly non-normal (anderson-darling test) 
        [pval_temp, h, ci_temp] = ranksum(temp_nosent, temp_sent, 'tail', 'right'); %Since they are highly non-normal (anderson-darling test)
        ci_temp = [0,0];

        pval_f = pval_f + pval_temp;

        if isempty(ci_f)
            ci_f = ci_temp;
        else
            ci_f = ci_f + ci_temp;
        end
    end
    pval_f = pval_f/1000;
    ci_f = ci_f/1000;



    [h,pval, ci] = ttest2(slopes_nosent, slopes_sent, 'tail', 'both');
end

partLevel = 0.5;
t_95_nosent = log(partLevel)./slopes_nosent_f/60;
t_95_sent = log(partLevel)./slopes_sent_f/60;

t_95_nosent_med = log(partLevel)/nosent_med/60;
t_95_sent_med = log(partLevel)/sent_med/60;
%t_95_sent_mean_lb = log(partLevel)/(nosent_mean - ci_f(1))/60;
%t_95_sent_mean_ub = log(partLevel)/(nosent_mean - ci_f(2))/60;

%ub2 = (nosent_mean - ci_f(1))/nosent_mean;
m2 = (sent_med)/nosent_med;
%lb2 = (nosent_mean - ci_f(2))/nosent_mean;



subplot(1,3,3);
cats = [zeros(size(t_95_nosent,1),1);ones(size(t_95_sent,1),1)];
cats = categorical(cats,[0,1], {'no airsentry', 'airsentry'});
boxplot([t_95_nosent; t_95_sent],cats);
ylabel('minutes');
title('time taken to reach 50%');
ylim([0,120]);
if (~isempty(slopes_sent))
    text(0,110,['p=',num2str(pval_f), ', med ratio=', num2str(m2), 't_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
    text(0,100,['n_{proc}=', num2str(nnz(usesAirSentry)), '/', num2str(nnz(~usesAirSentry)), ':  ', num2str(numel(lengths_sent)), '/', num2str(numel(lengths_nosent))]);
else
    text(0,110,['t_{nosent}=', num2str(t_95_nosent_med), 'mins, t_{sent}=', num2str(t_95_sent_med), 'mins']);
end
    
%histogram(t_95_nosent, 100);
%hold on;
%histogram(t_95_sent, 50);
%hold off;

saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_upper', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_interprocslopes_upper', label, '.png']));
%%

%%%
tempTable = removevars(upperGITable, {'StudyNumber', 'procedureMax_all', 'procedureMax_sum', 'procedureMax_sum_v', 'procedureTot_all', 'procedureTot_sum', 'procedureTot_sum_v', 'procedureMean_all', 'procedureMean_sum', 'procedureMean_sum_v','nearestEvent_all', 'nearestEvent_sum', 'preProcedureMax_all','preProcedureMax_sum','preProcedureMax_sum_v', 'nearestEvent_sum_v', 'tDiff_all', 'tDiff_sum', 'tDiff_sum_v', 'diameters'});
%tempForest = TreeBagger(1000, tempTable, upperGITable.procedureMax_sum, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');
tempForest = TreeBagger(1000, tempTable, upperGITable.procedureTot_sum./upperGITable.procedureDuration, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

imp = tempForest.OOBPermutedPredictorDeltaError;
%imp(imp < 0) = 0;
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
bar(imp);
title('Upper GI Variable importance');
ylabel('Predictor importance estimates');
xlabel('Predictors');
h = gca;
h.XTickLabel = tempForest.PredictorNames;
h.XTick = 1:size(tempForest.PredictorNames,2);
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';

subplot(1,2,2);
oobErrorBaggedEnsemble = oobError(tempForest);
plot(oobErrorBaggedEnsemble)
xlabel('Number of grown trees');
ylabel('Out-of-bag MSE');

saveas(gcf,fullfile(resultsFolder,['fullproc_forestugi_', label, '.fig']));
saveas(gcf,fullfile(resultsFolder,['fullproc_forestugi_', label, '.png']));

diary off;

%% Lower GI
if ~isempty(lowerGITable)
    tempTable = removevars(lowerGITable, {'StudyNumber', 'procedureMax_all', 'procedureMax_sum', 'procedureMax_sum_v', 'procedureTot_all', 'procedureTot_sum', 'procedureTot_sum_v', 'procedureMean_all', 'procedureMean_sum', 'procedureMean_sum_v','nearestEvent_all', 'nearestEvent_sum', 'preProcedureMax_all','preProcedureMax_sum','preProcedureMax_sum_v', 'nearestEvent_sum_v', 'tDiff_all', 'tDiff_sum', 'tDiff_sum_v', 'diameters'});
    %tempForest = TreeBagger(1000, tempTable, lowerGITable.procedureMax_sum, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');
    tempForest = TreeBagger(1000, tempTable, lowerGITable.procedureTot_sum./lowerGITable.procedureDuration, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature');%, 'NumPredictorsToSample', 'all');

    imp = tempForest.OOBPermutedPredictorDeltaError;
    %imp(imp < 0) = 0;
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    bar(imp);
    title('Lower GI Variable importance');
    ylabel('Predictor importance estimates');
    xlabel('Predictors');
    h = gca;
    h.XTickLabel = tempForest.PredictorNames;
    h.XTick = 1:size(tempForest.PredictorNames,2);
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = 'none';

    subplot(1,2,2);
    oobErrorBaggedEnsemble = oobError(tempForest);
    plot(oobErrorBaggedEnsemble)
    xlabel('Number of grown trees');
    ylabel('Out-of-bag MSE');

    saveas(gcf,fullfile(resultsFolder,['fullproc_forestlgi_', label, '.fig']));
    saveas(gcf,fullfile(resultsFolder,['fullproc_forestlgi_', label, '.png']));

    a = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = objFun(p, pd_h, pd_nh)
    temp1 = paramci(pd_h,'Alpha',p);
    temp2 = paramci(pd_nh,'Alpha', p);
    err = abs(temp1(1) - temp2(2));
end