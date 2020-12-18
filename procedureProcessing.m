%% Process individual procedures
% Author: George Gordon
% Date: 18/12/2020


% Clear any old data
clc;
clear variables;
close all;

dataDir = ['C:\Users\george\OneDrive - The University of Nottingham\SAVE\'];

fileList_raw = dir(dataDir);
    
filterFun = @(x) regexpi(x, '^[0-9]{8}');
temp = cellfun(filterFun, {fileList_raw.name}, 'UniformOutput', false); 
fileList_raw = fileList_raw(~cellfun(@isempty,temp));

nFiles = size(fileList_raw,1);

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
        
        P = test(4); % TODO ensure dates match
        
        [data, datatimes, eventTimes, eventNames, avSampleTime, bg_data, fg_data,data_v] = loadAnnotatedData(Y,M,D,P);
        
        if isnan(data) % no valid events detected
            continue;
        end
        
        % Reference data
        firstEventIdx = find(datatimes >= eventTimes(1), 1);
        preProcedure = data(:,1:firstEventIdx-1);
        preProcedure_av = mean(preProcedure,2,'omitnan')/avSampleTime;
        
        compFun = @(x) strcmpi(x,'procedure starts');
        procStart_temp = find(cellfun(compFun,table2cell(eventNames)));
        if isempty(procStart_temp)
            procStartIdx = 1;
        else
            procStart_temp = procStart_temp(1); %somtimes theres are 2 procedure ends!
            procStartIdx = find(datatimes >= eventTimes(procStart_temp), 1);
        end
        compFun = @(x) strcmpi(x,'procedure ends');
        procEnd_temp = find(cellfun(compFun,table2cell(eventNames)));
        if isempty(procEnd_temp)
            procEndIdx = size(data,2);
        else
            procEnd_temp = procEnd_temp(1); %somtimes theres are 2 procedure ends!
            procEndIdx = find(datatimes >= eventTimes(procEnd_temp), 1);
        end
        
        datatimes_proc = datatimes(procStartIdx:procEndIdx-1);
        procedure = data(:,procStartIdx:procEndIdx-1);
        procedure_v = data_v(:, procStartIdx:procEndIdx-1);
        [procedure_max, maxIdx] = max(procedure,[],2);
        [procedure_max_sum, maxIdx_sum] = max(nansum(procedure,1),[],2);
        [procedure_max_sum_v, maxIdx_sum_v] = max(nansum(procedure_v,1),[],2);
        
        refTime = datatimes_proc(1);
        nearestEventIdxes = knnsearch(seconds(eventTimes - refTime),seconds(datatimes_proc(maxIdx) - refTime));
        
        for k = 1:size(nearestEventIdxes,1)
            temp = table2cell(eventNames(nearestEventIdxes(k),1));
            nearestEventNames{k} = temp{1};
            tDiff(k) = datatimes_proc(maxIdx(k)) - eventTimes(nearestEventIdxes(k));
        end
        
        nearestEventIdxes_sum = knnsearch(seconds(eventTimes - refTime),seconds(datatimes_proc(maxIdx_sum) - refTime));
        temp = table2cell(eventNames(nearestEventIdxes_sum(1),1));
        nearestEventNames_sum = temp{1};
        tDiff_sum = datatimes_proc(maxIdx_sum) - eventTimes(nearestEventIdxes_sum);
            
        nearestEventIdxes_sum_v = knnsearch(seconds(eventTimes - refTime),seconds(datatimes_proc(maxIdx_sum_v) - refTime));
        temp = table2cell(eventNames(nearestEventIdxes_sum_v(1),1));
        nearestEventNames_sum_v = temp{1};
        tDiff_sum_v = datatimes_proc(maxIdx_sum_v) - eventTimes(nearestEventIdxes_sum_v);
        
        %% Fix should load variables to check
        compFun = @(x) strcmpi(x,'PR exam');
        prexam = find(cellfun(compFun,table2cell(eventNames)));
        if isempty(prexam)
            isUGI = false;
            nearestEventNames{:}
            nearestEventNames_sum
            nearestEventNames_sum_v
        else
            isUGI = true;
        end
        
        
    end
    
end