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
            procStartIdx = find(datatimes >= (eventTimes(procStart_temp) - minutes(0)), 1);
        end
        compFun = @(x) strcmpi(x,'procedure ends');
        procEnd_temp = find(cellfun(compFun,table2cell(eventNames)));
        if isempty(procEnd_temp)
            procEndIdx = size(data,2);
        else
            procEnd_temp = procEnd_temp(1); %somtimes theres are 2 procedure ends!
            procEndIdx = find(datatimes >= eventTimes(procEnd_temp), 1);
        end
        
%         compFun = @(x) strcmpi(x,'oesophageal extubation');
%         t = find(cellfun(compFun,table2cell(eventNames)));
%         if (~isempty(t))
%             a = 1;
%         end
        
        datatimes_proc = datatimes(procStartIdx:procEndIdx-1);
        procedure = data(:,procStartIdx:procEndIdx-1);
        procedure_v = data_v(:, procStartIdx:procEndIdx-1);
        [procedure_max, maxIdx] = max(procedure,[],2);
        [procedure_max_sum, maxIdx_sum] = max(nansum(procedure,1),[],2);
        [procedure_max_sum_v, maxIdx_sum_v] = max(nansum(procedure_v,1),[],2);
        
        procedureTot = sum(procedure,2);
        procedureTot_sum = sum(nansum(procedure,1),2);
        procedureTot_sum_v = sum(nansum(procedure_v,1),2);
        
        procedureMean = mean(procedure,2)/avSampleTime;
        procedureMean_sum = mean(nansum(procedure,1),2)/avSampleTime;
        procedureMean_sum_v = mean(nansum(procedure_v,1),2)/avSampleTime;
        
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
        
        tempResults.StudyNumber = otherVars.StudyNumber;
        
        tComb = join(otherVars,tempResults);
        
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
                
                count_5min_all = bg_next(:,idx_5min)/avSampleTime;
                count_5min_sum = sum(bg_next(:,idx_5min)/avSampleTime,1);

                count_10min_all = bg_next(:,idx_10min)/avSampleTime;
                count_10min_sum = sum(bg_next(:,idx_10min)/avSampleTime,1);

                count_20min_all = bg_next(:,idx_20min)/avSampleTime;
                count_20min_sum = sum(bg_next(:,idx_20min)/avSampleTime,1);

                tempResults_next.count_5min_all = {count_5min_all./(bg_next(:,1)/avSampleTime)};
                tempResults_next.count_5min_sum = count_5min_sum./(nansum(bg_next(:,1))/avSampleTime);
                tempResults_next.count_10min_all = {count_10min_all./(bg_next(:,1)/avSampleTime)};
                tempResults_next.count_10min_sum = count_10min_sum./(nansum(bg_next(:,1))/avSampleTime);
                tempResults_next.count_20min_all = {count_20min_all./(bg_next(:,1)/avSampleTime)};
                tempResults_next.count_20min_sum = count_20min_sum./(nansum(bg_next(:,1))/avSampleTime);
                
%                 tempResults_next.count_5min_all = {count_5min_all};
%                 tempResults_next.count_5min_sum = count_5min_sum;
%                 tempResults_next.count_10min_all = {count_10min_all};
%                 tempResults_next.count_10min_sum = count_10min_sum;
%                 tempResults_next.count_20min_all = {count_20min_all};
%                 tempResults_next.count_20min_sum = count_20min_sum;

                %TODO repeat for each 
                logvals = log(nansum(bg_next./repmat(bg_next(:,1),1,size(bg_next,2)),1));
                logdiff = diff(logvals)/avSampleTime;
                validslope = logdiff < 0;
                tempFiltBySlope = logvals;
                tempFiltBySlope(~validslope) = NaN;
                
                tempSlope = [];
                lengths = [];
                slopes = [];
                for kk = 1:size(tempFiltBySlope,2)
                    if isnan(tempFiltBySlope(kk))
                        if ~isempty(tempSlope)
                            x = 0:size(tempSlope,1)-1;
                            x = x*avSampleTime;
                            p = polyfit(x, tempSlope,1);
                            slopes = [slopes; p(1)];
                            lengths = [lengths; size(tempSlope,1)];
                            tempSlope = [];
                        end
                    else
                        tempSlope = [tempSlope; tempFiltBySlope(kk)];
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
                
                subplot(3,2,1);
                plot(minutes(opTime2_next - opTime2_next(1)), log(nansum(bg_next./bg_next(:,1),1)));
                xlim([0,30]);
                xlabel('Minutes after procedure end');
                hold on;
                
                subplot(3,2,3);

                plot(minutes(opTime2_next(1:end-1) - opTime2_next(1)), logdiff);
                xlim([0,30]);
                xlabel('Minutes after procedure end');
                hold on;
                
                subplot(3,2,5)
                plot(minutes(opTime2_next - opTime2_next(1)), tempFiltBySlope);
                xlim([0,30]);
                xlabel('Minutes after procedure end');
                hold on;
                
                subplot(3,2,2)
                histogram(slopes_temp,linspace(-0.01,0,30));
                xlabel('Slopes');
                pause(0.1);
                
                
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

maxPart = procedureTable.procedureMax_sum;
maxPart_cat = procedureTable.Procedure;

meanPart = procedureTable.procedureMean_sum;
meanPart_cat = procedureTable.Procedure;

figure;
subplot(1,3,1);
boxplot(totPart, totPart_cat);
title('Total no. particles');

subplot(1,3,2);
boxplot(maxPart, maxPart_cat);
title('Max no. particles/m^3s^{-1}');

subplot(1,3,3);
boxplot(meanPart, meanPart_cat);
title('Mean no. particles/m^3s^{-1}');

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
cats = [zeros(size(slopes_nosent_f,1),1);ones(size(slopes_sent_f,1),1)];
cats = categorical(cats,[0,1], {'airsentry', 'no airsentry'});
boxplot([slopes_nosent_f; slopes_sent_f], cats);

[h,pval] = ttest2(slopes_nosent_f, slopes_sent_f, 'tail', 'right');

partLevel = 0.01;
t_95_nosent = log(partLevel)./slopes_nosent_f/60;
t_95_sent = log(partLevel)./slopes_sent_f/60;

t_95_nosent_mean = log(partLevel)/nosent_mean/60;
t_95_sent_mean = log(partLevel)/sent_mean/60;

subplot(1,3,3);
histogram(t_95_nosent, 100);
hold on;
histogram(t_95_sent, 50);
hold off;

%%%
tempTable = removevars(upperGITable, {'StudyNumber', 'procedureMax_all', 'procedureMax_sum', 'procedureMax_sum_v', 'procedureTot_all', 'procedureTot_sum', 'procedureTot_sum_v', 'procedureMean_all', 'procedureMean_sum', 'procedureMean_sum_v','nearestEvent_all', 'nearestEvent_sum', 'nearestEvent_sum_v', 'tDiff_all', 'tDiff_sum', 'tDiff_sum_v'});
tempForest = TreeBagger(1000, tempTable, upperGITable.procedureMax_sum, 'CategoricalPredictors',[1,4:size(tempTable,2)], 'Method', 'Regression', 'OOBPrediction','On', 'OOBPredictorImportance','on','Surrogate','on');

imp = tempForest.OOBPermutedPredictorDeltaError;
%imp(imp < 0) = 0;
figure;
subplot(1,2,1);
bar(imp);
title('Variable importance');
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