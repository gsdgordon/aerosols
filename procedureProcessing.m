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

procedureTable = [];
interProcedureTable = [];

limitSize = false;
lt = true;
sizeLim = 5;

includeThroatSpray = true;
supressPosChanges = true;
fgOnly = true;

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
        
        [data, datatimes, eventTimes, eventNames, avSampleTime, otherVars, bg_data, fg_data, data_v, bg_data_v, fg_data_v, diameters, diameters_av, data_next, opTime2_next,  bg_data_next, fg_data_next, data_v_next, bg_data_v_next, fg_data_v_next] = loadAnnotatedData(Y,M,D,P);
        
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
    
        if ~isempty(regexpi(firstEvent, '.*doors locked.*'))
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

        % Supress position changes
        if supressPosChanges
            compFun = @(x) regexpi(x,'.*position changes');
            posChange_temp = cellfun(compFun,table2cell(eventNames),'UniformOutput', false);
            posChange_temp = cellfun(@isempty,posChange_temp);
            posChange_temp = ~posChange_temp;
            posChange_temp = find(posChange_temp);

            supressionWindow = 100;
            for posChange = posChange_temp.'
                startPosChIdx = find(datatimes >= (eventTimes(posChange)), 1);
                endPosChIdx = find(datatimes > (eventTimes(posChange) + seconds(supressionWindow)), 1) - 1;
                
                endRefIdx = startPosChIdx - 5;
                startRefIdx = find(datatimes > (eventTimes(posChange) - seconds(supressionWindow)), 1);
                
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
        
        if supressPosChanges
            tempResults.posChangesSupressed = posChangesSuppressed;
        end
        
        tComb = join(otherVars,tempResults);
        
        if (tComb.RoomType == 'theatre') || (tComb.PatientMask == 'yes')
            %Exclude masks and ERCP procedures
            continue;
        end
        
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

                %TODO repeat for each 
                logvals_norm = log(nansum(bg_next(validDiams,:),1)./nansum(bg_next(validDiams,1)));
                logvals = log(nansum(bg_next(validDiams,:),1)); %Do we need to divide by average sample time?
                
                logdiff = diff(logvals)/avSampleTime;
                logdiff_norm = diff(logvals_norm)/avSampleTime;
                validslope = logdiff < 0;
                tempFiltBySlope = logvals;
                tempFiltBySlope(~validslope) = NaN;
                
                tempFiltBySlope_norm = logvals_norm;
                tempFiltBySlope_norm(~validslope) = NaN;
                
                tempSlope = [];
                lengths = [];
                slopes = [];
                count = 0;
                sizethresh = 4; % Max length of segment to fit slope to
                for kk = 1:size(tempFiltBySlope,2)
                    if isnan(tempFiltBySlope(kk)) || count >= sizethresh
                        if ~isempty(tempSlope)
                            x = 0:size(tempSlope,1)-1;
                            x = x*avSampleTime;
                            p = polyfit(x, tempSlope,1);
                            slopes = [slopes; p(1)];
                            lengths = [lengths; size(tempSlope,1)];
                            tempSlope = [];
                        end
                        count = 0;
                    else
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
lowerGITable = procedureTable(procedureTable.Procedure == 'lower GI', :);

figure;
subplot(1,2,1);
pie(categorical(lowerGITable.nearestEvent_sum));
title('Lower GI max events');

subplot(1,2,2);
pie(categorical(upperGITable.nearestEvent_sum));
title('Upper GI max events');

figure;
subplot(1,2,1);
pie(categorical(lowerGITable.nearestEvent_sum_v));
title('Lower GI max events (volume)');

subplot(1,2,2);
pie(categorical(upperGITable.nearestEvent_sum_v));
title('Upper GI max events (volume)');

%%%
totPart = procedureTable.procedureTot_sum;
totPart_cat = procedureTable.Procedure;

totPart_pre = procedureTable.preProcedureAv_sum .* procedureTable.procedureDuration * 60;
totPart_pre_cat = categorical(ones(size(totPart_pre)),[0,1],{'a','pre-procedure'});

[~, p_lp, ci_lp] = ttest(log(totPart(totPart_cat == 'lower GI')), log(totPart_pre(totPart_cat == 'lower GI')));
[~, p_up, ci_up] = ttest(log(totPart(totPart_cat == 'upper GI')), log(totPart_pre(totPart_cat == 'upper GI')));
[~, p_lu, ci_lu] = ttest2(log(totPart(totPart_cat == 'lower GI')./procedureTable.procedureDuration(totPart_cat == 'lower GI')), log(totPart(totPart_cat == 'upper GI')./procedureTable.procedureDuration(totPart_cat == 'upper GI')));

ratio_lp = exp(ci_lp);
ratio_up = exp(ci_up);
ratio_lu = exp(ci_lu);

meanLGI = exp(mean(log(totPart(totPart_cat == 'lower GI'))));
meanUGI = exp(mean(log(totPart(totPart_cat == 'upper GI'))));

mean_lp = sqrt(ratio_lp(1) * ratio_lp(2));
mean_up = sqrt(ratio_up(1) * ratio_up(2));
mean_lu = sqrt(ratio_lu(1) * ratio_lu(2));

disp(['Total LGI-pre ratio: ', num2str(mean_lp), ' (', num2str(ratio_lp(1)), '-', num2str(ratio_lp(2)), ') p=', num2str(p_lp), ', mean = ', num2str(meanLGI)]);
disp(['Total UGI-pre ratio: ', num2str(mean_up), ' (', num2str(ratio_up(1)), '-', num2str(ratio_up(2)), ') p=', num2str(p_up), ', mean = ', num2str(meanUGI)]);
disp(['Total LGI-UGI ratio: ', num2str(mean_lu), ' (', num2str(ratio_lu(1)), '-', num2str(ratio_lu(2)), ') p=', num2str(p_lu)]);

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

[~, p_lp, ci_lp] = ttest2(log(maxPartPre), log(maxPart(maxPart_cat == 'lower GI')));
[~, p_up, ci_up] = ttest2(log(maxPartPre), log(maxPart(maxPart_cat == 'upper GI')));
[~, p_lu, ci_lu] = ttest2(log(maxPart(maxPart_cat == 'lower GI')), log(maxPart(maxPart_cat == 'upper GI')));

ratio_lp = exp(ci_lp);
ratio_up = exp(ci_up);
ratio_lu = exp(ci_lu);

mean_lp = sqrt(ratio_lp(1) * ratio_lp(2));
mean_up = sqrt(ratio_up(1) * ratio_up(2));
mean_lu = sqrt(ratio_lu(1) * ratio_lu(2));

meanPart = procedureTable.procedureMean_sum;
meanPart_cat = procedureTable.Procedure;

figure;
subplot(1,3,1);
%boxplot(totPart, totPart_cat);
violinplot([totPart./totPart_pre], [totPart_cat]);
title('Total no. particles');

subplot(1,3,2);
%boxplot([maxPartPre; maxPart; maxPart_diff], [maxPartPre_cat; maxPart_cat; maxPart_diffCat]);
violinplot([maxPartPre; maxPart], [maxPartPre_cat; maxPart_cat]);
title('Max no. particles during proc./m^3s^{-1}');

subplot(1,3,3);
%boxplot(meanPart, meanPart_cat);
violinplot(meanPart, meanPart_cat);
title('Mean no. particles/m^3s^{-1}');

%% Durations
durations = procedureTable.procedureDuration;
durations_LGI = durations(procedureTable.Procedure == 'lower GI');
durations_UGI = durations(procedureTable.Procedure == 'upper GI');

figure;
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

disp(['UGI proc. duration: Mean: ', num2str(mean_duration_UGI), ', Med: ', num2str(med_duration_UGI), ' mins']);
disp(['LGI proc. duration: Mean: ', num2str(mean_duration_LGI), ', Med: ', num2str(med_duration_LGI), ' mins']);


figure;
subplot(1,2,1);
scatter(lowerGITable.procedureDuration,lowerGITable.procedureTot_sum);
[fitresult, gof, ~] = fit(lowerGITable.procedureDuration,lowerGITable.procedureTot_sum,'poly1');
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

subplot(1,2,2);
scatter(upperGITable.procedureDuration,upperGITable.procedureTot_sum);
[fitresult, gof, ~] = fit(upperGITable.procedureDuration,upperGITable.procedureTot_sum,'poly1');
newx = linspace(min(upperGITable.procedureDuration), max(upperGITable.procedureDuration),100);
yfit = feval(fitresult,newx);

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);

if (size(upperGITable.procedureTot_sum,1) > 2)
    p21 = predint(fitresult,newx,0.95,'functional','off');
end
hold on;
plot(newx, yfit, 'k');

if (size(upperGITable.procedureTot_sum,1) > 2)
    plot(newx, p21, 'm--');
end
text(min(newx)*1.1, max(upperGITable.procedureTot_sum)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);
hold off;
xlabel('procecure length (mins)');
ylabel('tot no. particles');
title('upper GI tot. particles v length');



%% Plot vars

%% Anal tone
figure;
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

text(1.0, max(maxPart)*0.8,['Low-Medium: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);
text(1.0, max(maxPart)*0.6,['Low-High: ', num2str(sqrt(ci_lh(1)*ci_lh(2))), ' (', num2str(ci_lh(1)), '-', num2str(ci_lh(2)), ') p = ', num2str(p_lh)]);
text(1.0, max(maxPart)*0.4,['Medium-High: ', num2str(sqrt(ci_mh(1)*ci_mh(2))), ' (', num2str(ci_mh(1)), '-', num2str(ci_mh(2)), ') p = ', num2str(p_mh)]);

%% Sedation, proc. duration
figure;
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

figure;
sedation = procedureTable.Sedation;
sedDuration = procedureTable.procedureDuration;

sedation = sedation(procedureTable.Procedure == 'lower GI');
sedDuration = sedDuration(procedureTable.Procedure == 'lower GI');
violinplot(sedDuration, sedation);
title('Duration particles vs. sedation, LGI');
ylabel('duration (mins)');

[~, p_sed, ci_sed] = ttest2(log(sedDuration(sedation == 'midazolam')), log(sedDuration(sedation == 'entonox')));
ci_sed = exp(ci_sed);

text(1.0, max(sedDuration)*0.8,['midazolam-entonox: ', num2str(sqrt(ci_sed(1)*ci_sed(2))), ' (', num2str(ci_sed(1)), '-', num2str(ci_sed(2)), ') p = ', num2str(p_sed)]);

%%
figure;
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

figure;
sex = procedureTable.Sex;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sex = sex(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, sex);
title('Max particles vs. sex lower GI');
ylabel('max no particles /m^3/s');

[~, p_mf, ci_mf] = ttest2(log(maxPart(sex == 'male')), log(maxPart(sex == 'female')));
ci_mf = 1./exp(ci_mf);

text(1.0, max(maxPart)*0.8,['Male-female: ', num2str(sqrt(ci_mf(1)*ci_mf(2))), ' (', num2str(ci_mf(1)), '-', num2str(ci_mf(2)), ') p = ', num2str(p_mf)]);

%% Sedation
figure;
sedation = procedureTable.Sedation;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sedation = sedation(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, sedation);
title('Max particles vs. sed upper');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(sedation == 'midazolam')), log(maxPart(sedation == 'none')));
ci_mn = exp(ci_mn);

text(1.0, max(maxPart)*0.8,['Midaz-none: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);

figure;
sedation = procedureTable.Sedation;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

sedation = sedation(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, sedation);
title('Max particles vs. sed lower');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(sedation == 'midazolam')), log(maxPart(sedation == 'entonox')));

ci_mn = exp(ci_mn);

text(1.0, max(maxPart)*0.8,['Midaz-entonox: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);

%% Smoker
figure;
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


figure;
smoker = procedureTable.Smoker;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

smoker = smoker(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, smoker);
title('Max particles vs. smoker');
ylabel('max no particles /m^3/s');

[~, p_mn, ci_mn] = ttest2(log(maxPart(smoker == 'yes')), log(maxPart(smoker == 'no')));

ci_mn = exp(ci_mn);
text(1.0, max(maxPart)*0.8,['Smoker-nonsmoker: ', num2str(sqrt(ci_mn(1)*ci_mn(2))), ' (', num2str(ci_mn(1)), '-', num2str(ci_mn(2)), ') p = ', num2str(p_mn)]);

%% Discomfort
figure;
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

figure;
discomfort = procedureTable.Discomfort;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

discomfort = discomfort(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
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


%% Co2
figure;
CO2 = procedureTable.UseOfCO2orWater;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

CO2 = CO2(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, CO2);
title('Max particles vs. use of CO2');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(CO2 == 'CO2')), log(maxPart(CO2 == 'Water')));

ci_lm = exp(ci_lm);
text(1.0, max(maxPart)*0.8,['CO2-water: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);


%% Hysterectomy
figure;
hyst = procedureTable.PreviousHysterectomy;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

hyst = hyst(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');
violinplot(maxPart, hyst);
title('Max particles vs. previous hyst');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(hyst == 'yes')), log(maxPart(hyst == 'no')));

ci_lm = 1./exp(ci_lm);
text(1.0, max(maxPart)*0.8,['Yes-No: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);


%% UGI route
figure;
ugiroute = procedureTable.UGIroute;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

ugiroute = ugiroute(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, ugiroute);
title('Max particles vs. UGI route');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(ugiroute == 'nasal')), log(maxPart(ugiroute == 'oral')));

ci_lm = exp(ci_lm);
text(1.0, max(maxPart)*0.8,['Nasal-Oral: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);


%% Hernia
figure;
hiatus = procedureTable.HiatusHernia;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

hiatus = hiatus(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, hiatus);
title('Max particles vs. hiatus');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(hiatus == 'yes')), log(maxPart(hiatus == 'unknown')));

ci_lm = exp(ci_lm);
text(1.0, max(maxPart)*0.8,['Yes-No: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);


%% Suctioning
figure;
suctioning = procedureTable.Suctioning;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;

suctioning = suctioning(procedureTable.Procedure == 'upper GI');
maxPart = maxPart(procedureTable.Procedure == 'upper GI');
violinplot(maxPart, suctioning);
title('Max particles vs. suctioning');
ylabel('max no particles /m^3/s');

[~, p_lm, ci_lm] = ttest2(log(maxPart(suctioning == 'yes')), log(maxPart(suctioning == 'unknown')));

ci_lm = exp(ci_lm);
text(1.0, max(maxPart)*0.8,['Yes-no: ', num2str(sqrt(ci_lm(1)*ci_lm(2))), ' (', num2str(ci_lm(1)), '-', num2str(ci_lm(2)), ') p = ', num2str(p_lm)]);



%% Age
figure;
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
text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);


figure;
age = procedureTable.Age;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
age = age(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');

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
a = 1;

confidenceInt  = confint(fitresult);
slopeconf = confidenceInt(:,1);
text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);



%% BMI
figure;
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
text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);


figure;
BMI = procedureTable.BMI;
maxPart = procedureTable.procedureTot_sum ./ procedureTable.procedureDuration * 20;
BMI = BMI(procedureTable.Procedure == 'lower GI');
maxPart = maxPart(procedureTable.Procedure == 'lower GI');

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
text(min(newx)*1.1, max(maxPart)*0.8,['y = ', num2str(fitresult.p1), '(', num2str(min(slopeconf)), ' - ', num2str(max(slopeconf)), ')x + ', num2str(fitresult.p2), ', r = ', num2str(gof.rsquare)]);



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

figure;
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

% Now compute slopes
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

figure;
subplot(1,3,1)
histogram(slopes_nosent_f,50);
nosent_mean = sum(lengths_nosent/sum(lengths_nosent) .* slopes_nosent);
sent_mean = sum(lengths_sent/sum(lengths_sent) .* slopes_sent);
xlabel('slope');
ylabel('frequency');
title(['slopes w/o and w airsentry (means = ', num2str(nosent_mean), ', ', num2str(sent_mean), ')']);

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

pval_f = 0;
for k=1:1000
    temp_nosent = datasample(slopes_nosent,size(slopes_nosent,1),'Weights',lengths_nosent);
    temp_sent = datasample(slopes_sent,size(slopes_sent,1)*2,'Weights',lengths_sent);

    [h,pval_temp, ci_temp] = ttest2(temp_nosent, temp_sent, 'tail', 'both');
    
    pval_f = pval_f + pval_temp;
end
pval_f = pval_f/1000;

[h,pval, ci] = ttest2(slopes_nosent, slopes_sent, 'tail', 'both');

partLevel = 0.42;
t_95_nosent = log(partLevel)./slopes_nosent_f/60;
t_95_sent = log(partLevel)./slopes_sent_f/60;

t_95_nosent_mean = log(partLevel)/nosent_mean/60;
t_95_sent_mean = log(partLevel)/sent_mean/60;
t_95_sent_mean_lb = log(partLevel)/(nosent_mean - ci(1))/60;
t_95_sent_mean_ub = log(partLevel)/(nosent_mean - ci(2))/60;

subplot(1,3,3);
histogram(t_95_nosent, 100);
hold on;
histogram(t_95_sent, 50);
hold off;

%%%
tempTable = removevars(upperGITable, {'StudyNumber', 'procedureMax_all', 'procedureMax_sum', 'procedureMax_sum_v', 'procedureTot_all', 'procedureTot_sum', 'procedureTot_sum_v', 'procedureMean_all', 'procedureMean_sum', 'procedureMean_sum_v','nearestEvent_all', 'nearestEvent_sum', 'preProcedureMax_all','preProcedureMax_sum','preProcedureMax_sum_v', 'nearestEvent_sum_v', 'tDiff_all', 'tDiff_sum', 'tDiff_sum_v', 'diameters'});
%tempForest = TreeBagger(1000, tempTable, upperGITable.procedureMax_sum, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');
tempForest = TreeBagger(1000, tempTable, upperGITable.procedureTot_sum./upperGITable.procedureDuration, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');

imp = tempForest.OOBPermutedPredictorDeltaError;
%imp(imp < 0) = 0;
figure;
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

%% Lower GI
tempTable = removevars(lowerGITable, {'StudyNumber', 'procedureMax_all', 'procedureMax_sum', 'procedureMax_sum_v', 'procedureTot_all', 'procedureTot_sum', 'procedureTot_sum_v', 'procedureMean_all', 'procedureMean_sum', 'procedureMean_sum_v','nearestEvent_all', 'nearestEvent_sum', 'preProcedureMax_all','preProcedureMax_sum','preProcedureMax_sum_v', 'nearestEvent_sum_v', 'tDiff_all', 'tDiff_sum', 'tDiff_sum_v', 'diameters'});
%tempForest = TreeBagger(1000, tempTable, lowerGITable.procedureMax_sum, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature', 'NumPredictorsToSample', 'all');
tempForest = TreeBagger(1000, tempTable, lowerGITable.procedureTot_sum./lowerGITable.procedureDuration, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on','PredictorSelection', 'interaction-curvature');%, 'NumPredictorsToSample', 'all');

imp = tempForest.OOBPermutedPredictorDeltaError;
%imp(imp < 0) = 0;
figure;
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

a = 1;