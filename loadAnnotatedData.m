function [data, opTime2, eventTimes, eventNames, avSampleTime, T_othervars, bg_data, fg_data, data_v, bg_data_v, fg_data_v, diameters, diameters_av, data_next, opTime2_next,  bg_data_next, fg_data_next, data_v_next, bg_data_v_next, fg_data_v_next] = loadAnnotatedData(Y,M,D,P)

    datestring = sprintf('%0.4d%0.2d%0.2d',Y,M,D);

    folder = ['C:\Users\george\OneDrive - The University of Nottingham\SAVE\',datestring,'\'];

    if exist(fullfile(folder, [datestring,'_aerotrak.xlsx']))
        T_raw = readtable(fullfile(folder, [datestring,'_aerotrak.xlsx']));
        %T = readtable('/home/george/Desktop/20201030_aerotrak.xlsx');

        T_procindex_raw = readtable('C:\Users\george\OneDrive - The University of Nottingham\SAVE\ProcedureIndex.xlsx', 'ReadVariableNames', true, 'Format', 'auto');

        annotationFile = fullfile(folder,[datestring, '_patient', num2str(P), '.csv']);
        T_annotation = readtable(annotationFile, 'ReadVariableNames', false, 'HeaderLines', 0, 'Format', 'auto'); % FIX ignores first row

        %% Exclude if necessary
        exclude = false;
        temp = T_annotation(strcmpi(T_annotation.Var1, 'Study number'),2);
        temp = table2cell(temp);
        temp = temp{1};

        if isnan(temp)
            exclude = true;
        else
            tempStudyNo = str2num(temp);
            
            tableIdx = find(T_procindex_raw.StudyNo_ == tempStudyNo);

            exclude_temp = T_procindex_raw.Exclude(tableIdx);
            exclude_temp = exclude_temp{1};

            temp = regexp(exclude_temp, 'y.*');

            if ~isempty(temp)
                exclude = true;
            end
        end
    else
        exclude = true;
    end
    
    if exclude
        % Exclude this data point
        data = NaN;
        opTime2 = NaN;
        eventTimes = NaN;
        eventNames = NaN;
        avSampleTime = NaN;
        bg_data = NaN;
        fg_data = NaN;
        data_v = NaN;
        bg_data_v = NaN;
        fg_data_v = NaN;
        diameters = NaN;
        diameters_av = NaN;
        
        T_othervars = NaN;
        
        data_next = NaN;
        opTime2_next = NaN;
        bg_data_next = NaN;
        fg_data_next = NaN;
        data_v_next = NaN;
        bg_data_v_next = NaN;
        fg_data_v_next = NaN;
        return;
    end
    
    %%

    opTime = T_raw.DateAndTime;
    location = T_raw.Location;
    sampleTime = T_raw.SampleTime;
    avSampleTime = mode(sampleTime); % assumes sample time is not changed during operation, but excludes partial samples


    
    nextAnnotationFile = fullfile(folder,[datestring, '_patient', num2str(P+1), '.csv']);
    
    nextPatientExists = false;
    if exist(nextAnnotationFile)
        nextPatientExists = true;
        T_nextannotation = readtable(nextAnnotationFile, 'ReadVariableNames', false, 'HeaderLines', 0, 'Format', 'auto'); % FIX ignores first row
    end

    foundEventStart = false;
    validEvents = false(size(T_annotation,1),1);
    for k=1:size(T_annotation,1)
        currentText = table2array(T_annotation(k,1));
        currentText = currentText{1};
        if strcmpi(currentText, 'incident')
            eventStart_idx = k+1;
            foundEventStart = true;
            continue
        end

        currentTime = table2array(T_annotation(k,2));
        if foundEventStart
            if ~isempty(currentText)
                if ~(currentText(end) == '*' && currentText(end-1) == "*")
                    validEvents(k) = true;
                end
            end
        end
    end
    
    if nnz(validEvents) == 0
        data = NaN;
        opTime2 = NaN;
        eventTimes = NaN;
        eventNames = NaN;
        avSampleTime = NaN;
        bg_data = NaN;
        fg_data = NaN;
        data_v = NaN;
        bg_data_v = NaN;
        fg_data_v = NaN;
        diameters = NaN;
        diameters_av = NaN;
        
        T_othervars = NaN;
        
        data_next = NaN;
        opTime2_next = NaN;
        bg_data_next = NaN;
        fg_data_next = NaN;
        data_v_next = NaN;
        bg_data_v_next = NaN;
        fg_data_v_next = NaN;
        return;
    end

    eventNames = T_annotation(validEvents,1);
    %eventTimes_endo = T_annotation(validEvents,2);
    %eventTimes_obscam = T_annotation(validEvents,3);
    temp_times = table2array(T_annotation(validEvents,4));
    eventTimes_aerotrak = datetime(temp_times);
    
    if nextPatientExists
        foundEventStart_next = false;
        validEvents_next = false(size(T_nextannotation,1),1);
        for k=1:size(T_nextannotation,1)
            currentText = table2array(T_nextannotation(k,1));
            currentText = currentText{1};
            if strcmpi(currentText, 'incident')
                eventStart_idx_next = k+1;
                foundEventStart_next = true;
                continue
            end

            currentTime = table2array(T_nextannotation(k,2));
            if foundEventStart_next
                if ~isempty(currentText)
                    if ~(currentText(end) == '*' && currentText(end-1) == "*")
                        validEvents_next(k) = true;
                    end
                end
            end
        end
        
        eventNames_next = T_nextannotation(eventStart_idx_next:end,1);
        temp_times = table2array(T_nextannotation(eventStart_idx_next:end,4));
        eventTimes_aerotrak_next = datetime(temp_times);
    
        nLoops = 2;
    else
        nLoops = 1;
    end
    
    %% Now load and process other variables
    T_othervars_raw = readtable(annotationFile, 'ReadVariableNames', false, 'HeaderLines', 0, 'Format', 'auto');
    T_othervars_raw = T_othervars_raw(1:eventStart_idx-2,1:2);
    T_othervars = table;
    
    nOtherVars = size(T_othervars_raw,1);
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Study number'),2);
    temp = table2cell(temp);
    T_othervars.StudyNumber = str2num(temp{1});
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'procedure'),2);
    temp = table2cell(temp);
    [procedure, ~] = processVars(temp, {'upper GI','lower GI', 'N/A'}, {'upper gi.*'}, {'lower gi.*'}, {''});
    T_othervars.Procedure = procedure;
    
        
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Age'),2);
    temp = table2cell(temp);
    T_othervars.Age = str2num(temp{1});
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'BMI'),2);
    temp = table2cell(temp);
    if isempty(temp{1})
        T_othervars.BMI = NaN;
    else
        T_othervars.BMI = str2num(temp{1});
    end
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'sedation used'),2);
    temp = table2cell(temp);
    [sedation, ~] = processVars(temp, {'midazolam','propofol', 'entonox', 'none'}, {'.*mid.*'}, {'.*prop.*'}, {'.*ento.*'},{'.*none.*','.*N/A.*',''});
    T_othervars.Sedation = sedation;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'sex'),2);
    temp = table2cell(temp);
    [sex, ~] = processVars(temp, {'male','female', 'unknown'}, {'male.*', 'm'}, {'female.*','f'},{'.*none.*','.*N/A.*',''});
    T_othervars.Sex = sex;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Anal tone'),2);
    temp = table2cell(temp);
    [analTone, ~] = processVars(temp, {'low','medium', 'high', 'N/A'}, {'.*low.*', '.*loose.*', '.*poor.*'}, {'.*normal.*', '.*medium.*', '.*middle.*'}, {'.*good.*', '.*tight.*', '.*high.*'},{'.*none.*','.*N/A.*',''});
    T_othervars.AnalTone = analTone;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Use of CO2 or water'),2);
    temp = table2cell(temp);
    [useOfCO2orWater, ~] = processVars(temp, {'CO2','Water', 'N/A'}, {'.*co2.*', '.*yes.*'}, {'.*water.*', '.*h2o.*'}, {'.*none.*','.*N/A.*',''});
    T_othervars.UseOfCO2orWater = useOfCO2orWater;
         
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Smoker'),2);
    temp = table2cell(temp);
    [isSmoker, ~] = processVars(temp, {'yes','no', 'unknown'}, {'.*occasionally.*', '.*yes.*', '.*[1-9]?[0-9]*.*'}, {'.*no.*', '.*[0]+.*'}, {'.*none.*','.*N/A.*',''});
    T_othervars.Smoker = isSmoker;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'previous hysterectomy'),2);
    temp = table2cell(temp);
    [previousHyst, ~] = processVars(temp, {'yes','no', 'unknown'}, {'y.*'}, {'.*no', 'n'}, {'.*none.*','.*N/A.*',''});
    T_othervars.PreviousHysterectomy = previousHyst;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'extensive diverticular disease'),2);
    temp = table2cell(temp);
    [diverticularDisease, ~] = processVars(temp, {'none', 'mild', 'extensive', 'unknown'}, {'.*none.*', '.*no.*'}, {'.*mild.*', '.*moderate.*'},{'.*extensive.*', '.*yes.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});    
    T_othervars.DiverticularDisease = diverticularDisease;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Route of UGI'),2);
    temp = table2cell(temp);
    [ugiRoute, ~] = processVars(temp, {'oral', 'nasal', 'N/A'}, {'.*oral.*', '.*mouth.*'}, {'.*nasal.*', '.*nose.*'}, {'.*none.*','.*N/A.*',''});
    ugiRoute = addcats(ugiRoute,'nasal abandonded');
    % Check if chase from nasal to oral
    ugiRoute_temp = T_procindex_raw.UGIRoute(tableIdx);
    temp = regexpi(ugiRoute_temp, '.*half.*');
    if ~isempty(temp{1})
        if ugiRoute == 'nasal'
            ugiRoute = 'nasal abandoned';
        end
    end
    T_othervars.UGIroute = ugiRoute;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Degree of looping'),2);
    temp = table2cell(temp);
    [looping, ~] = processVars(temp, {'low', 'medium', 'high', 'unknown/N/A'}, {'.*none.*', '.*low.*', '.*no.*'}, {'.*mild.*', '.*moderate.*', '.*medium.*'},{'.*extensive.*', '.*high.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});
    T_othervars.Looping = looping;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Patient discomfort'),2);
    temp = table2cell(temp);
    [discomfort, ~] = processVars(temp, {'low', 'medium', 'high', 'unknown'}, {'.*none.*', '.*low.*', '.*no.*'}, {'.*mild.*', '.*moderate.*', '.*medium.*'},{'.*extensive.*', '.*high.*', '.*severe.*'}, {'.*unknown.*','.*N/A.*',''});
    T_othervars.Discomfort = discomfort;
    
    % Hiatus hernia
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Hiatus hernia length(cm)'),2);
    temp = table2cell(temp);
    temp2 = temp{1};
    if (isnumeric(temp))
        hiatusHernia = categorical(~isnan(temp2)+1,1:2,{'no', 'yes'});
    else
        [hiatusHernia, ~] = processVars(temp, {'no', 'yes', 'unknown'}, {'.*none.*', '.*[0].*', '.*no.*'}, {'.*[1-9]?[0-9]*.*', '.*yes.*', '.*massive.*'}, {'.*unknown.*','.*N/A.*',''});
    end
    T_othervars.HiatusHernia = hiatusHernia;
    
    temp = T_othervars_raw(strcmpi(T_othervars_raw.Var1, 'Use of intermittent suctioning'),2);
    temp = table2cell(temp);
    [suctioning, ~] = processVars(temp, {'no', 'yes', 'unknown/NA'}, {'.*none.*', '.*no.*'}, {'.*yes.*', '.*y.*'}, {'.*unknown.*','.*N/A.*',''});
    T_othervars.Suctioning = suctioning;
    
    % Last 2 vars pull from proc. index
    %T_procindex_raw = readtable('C:\Users\george\OneDrive - The University of Nottingham\SAVE\ProcedureIndex.xlsx', 'ReadVariableNames', true, 'Format', 'auto');
    
    tableIdx = find(T_procindex_raw.StudyNo_ == T_othervars.StudyNumber);
    
    roomType_temp = T_procindex_raw.RoomType(tableIdx);
    [roomType, ~] = processVars(roomType_temp, {'endoscopy', 'theatre', 'unknown/NA'}, {'.*endoscopy.*'}, {'.*theatre.*'}, {'.*unknown.*','.*N/A.*',''});
    T_othervars.RoomType = roomType;
    
    patientMask_temp = T_procindex_raw.PatientMaskUsed(tableIdx);
    [patientMask, ~] = processVars(patientMask_temp, {'yes', 'no'}, {'.*yes.*', 'y'}, {'.*no.*', 'n','.*unknown.*','.*N/A.*',''});
    T_othervars.PatientMask = patientMask;
    
    airSentryUsed_temp = roomType_temp{1};
    temp = regexp(airSentryUsed_temp, '.*[Ss]+entry.*');
    if isempty(temp)
        T_othervars.AirSentryUsed = false;
    else
        T_othervars.AirSentryUsed = true;
    end
    
    %T_othervars = T_othervars';
    %T_othervars.RoomType = roomType;
    
    
%         % Mask
%         usesPatientMask_raw = ismember(patientNos,[59,60,61, 64,65]); %FIX should pull this
%         usesPatientMask = categorical(usesPatientMask_raw+1,1:2,{'No Mask', 'Mask'});
%         
%         % Room type
%         roomType_raw = ismember(patientNos,[31,32,33]); %FIX should pull this
%         roomType = categorical(roomType_raw+1,1:2,{'Endoscopy room', 'Theatre'});
%         
%           
    
    %%
    for k=1:size(eventNames,1)
        currentTime_dur =(eventTimes_aerotrak(k,1)) - datetime(year(eventTimes_aerotrak(k,1)), month(eventTimes_aerotrak(k,1)), day(eventTimes_aerotrak(k,1)),0,0,0); % FIX need to sort by date
        eventTimes(k,1) = datetime(Y,M,D,0,0,0) + currentTime_dur - seconds(avSampleTime); % Shift events, not aerotrak
    end
    
    [eventTimes, sortOrder] = sort(eventTimes,'ascend');
    eventNames = eventNames(sortOrder,1);
    
    if nextPatientExists
        for k=1:size(eventNames_next,1)
            currentTime_dur_next =(eventTimes_aerotrak_next(k,1))  - datetime(year(eventTimes_aerotrak_next(k,1)), month(eventTimes_aerotrak_next(k,1)), day(eventTimes_aerotrak_next(k,1)),0,0,0);; % FIX need to sort by date
            eventTimes_next(k,1) = datetime(Y,M,D,0,0,0) + currentTime_dur_next - seconds(avSampleTime); % Shift events, not aerotrak
        end
    
        [eventTimes_next, sortOrder_next] = sort(eventTimes_next,'ascend');
        eventNames_next = eventNames_next(sortOrder_next,1);
    end

    
    for n=1:nLoops
        
        if n == 1
            startTime = eventTimes(1,1) - duration(0,5,0);
            endTime = eventTimes(end,1) + duration(0,5,0);
        else % Between procedures
            compFun = @(x) strcmpi(x,'procedure ends');
            procEnd_temp = find(cellfun(compFun,table2cell(eventNames)));
            if isempty(procEnd_temp)
                startTime = eventTimes(end,1);
            else
                procEnd_temp = procEnd_temp(end); %somtimes theres are 2 procedure ends!
                startTime = eventTimes(procEnd_temp);
            end
            
            startTime = startTime + minutes(1);
            
            compFun = @(x) strcmpi(x,'procedure starts');
            procStart_temp = find(cellfun(compFun,table2cell(eventNames_next)));
            if isempty(procStart_temp)
                endTime = eventTimes_next(1,1);
            else
                procStart_temp = procStart_temp(1); %somtimes theres are 2 procedure ends!
                endTime = eventTimes_next(procStart_temp);
            end
            
            endTime = endTime - minutes(1);
            
            a = 1;
        end

        aeroOffsetTime = 0;

        tValid = isbetween(opTime,startTime,endTime);
        tValid = tValid & strcmpi(location,'Location01');

        T = T_raw(tValid,:);
        T = T(1:end-1,:); % remove last count as it is likely partial

        if n==1
            opTime2 = T.DateAndTime - aeroOffsetTime;
            opTime2_next = [];
        else
            if nnz(tValid) == 0
                continue;
            end
            opTime2_next = T.DateAndTime - aeroOffsetTime;
        end
        airVol = T.Volume_L_;

        % Now get the particle 'diameters' representing the edges of the counting
        % bins
        maxDiameter = 25; % in microns - this is from the spec sheet
        diameters = [T.Ch1Size__m_(1),T.Ch2Size__m_(1),T.Ch3Size__m_(1),T.Ch4Size__m_(1),T.Ch5Size__m_(1),T.Ch6Size__m_(1), maxDiameter];
        diameters_av = (diameters(1:end-1)+diameters(2:end))/2; % Mean diameter in each bin

        % The actual counts in each bin
        data_temp = [T.Ch1Diff___,T.Ch2Diff___,T.Ch3Diff___,T.Ch4Diff___,T.Ch5Diff___,T.Ch6Diff___];
        data_temp = data_temp./repmat(airVol,1,size(data_temp,2)) * 1000; %Because volume is liters so times by 1000 to get to m^3
        data_temp = data_temp';

        [bg_data_temp, fg_data_temp] = splitBGFG(data_temp, avSampleTime, true(1,size(data_temp,2)));

        % Convert counts to volumes
        vols = 4/3*pi*(diameters_av/2).^3 * (1e-6)^3;

        bin_sizes = diameters(2:7) - diameters(1:6);
        log_bin_sizes = log(diameters(2:7)) - log(diameters(1:6));

        data_v_temp = data_temp .* repmat(vols',1,size(data_temp,2));
        bg_data_v_temp = bg_data_temp .* repmat(vols',1,size(data_temp,2));
        fg_data_v_temp = fg_data_temp .* repmat(vols',1,size(data_temp,2));

        % Densities so that a probability density approach can be used
        data_v_density = data_v_temp ./ repmat(log_bin_sizes',1,size(data_temp,2)); %Try using log binsizes
        data_density = data_temp ./ repmat(log_bin_sizes',1,size(data_temp,2));
        %data_v_density = data_v ./ repmat(bin_sizes,size(data,1),1); %Try using linear binsizes
        %data_density = data ./ repmat(bin_sizes,size(data,1),1);

        bg_density = bg_data_temp ./ repmat(log_bin_sizes',1,size(data_temp,2));
        fg_density = fg_data_temp ./ repmat(log_bin_sizes',1,size(data_temp,2));

        nSizes = size(data_temp,1);
        tColor = lines(nSizes);

        %%
        useTubeCorrection = true;

        if useTubeCorrection
            tubeCorrection_tab = readtable('C:\Users\george\OneDrive - The University of Nottingham\SAVE\TubeCalibration\TubeBendCorrection.csv');
            %tubeCorrection_tab = readtable('/home/george/Desktop/TubeBendCorrection.csv');
            tubeCorrection = table2array(tubeCorrection_tab);

            for k=1:nSizes
                correctionIdx = find(diameters(k) == tubeCorrection(:,1));

                correctionVal(k) = tubeCorrection(correctionIdx,2);
            end
        else
            correctionVal = zeros(1,nSizes);
        end

        data_temp = data_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        bg_data_temp = bg_data_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        fg_data_temp = fg_data_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        data_v_temp = data_v_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        bg_data_v_temp = bg_data_v_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        fg_data_v_temp = fg_data_v_temp ./ repmat(1-correctionVal',1,size(data_temp,2));
        
        if n == 1
            data = data_temp;
            bg_data = bg_data_temp;
            fg_data = fg_data_temp;
            data_v = data_v_temp;
            bg_data_v = bg_data_v_temp;
            fg_data_v = fg_data_v_temp;
            
            data_next = [];
            bg_data_next = [];
            fg_data_next = [];
            data_v_next = [];
            bg_data_v_next = [];
            fg_data_v_next = [];
            opTime2_next = [];
        else
            tDiff = diff(opTime2_next);
            maxTdiff = max(tDiff);
            if maxTdiff > 2*seconds(avSampleTime)
                validInterProc = false;
            else
                validInterProc = true;
            end
            
            if validInterProc
                data_next = data_temp;
                bg_data_next = bg_data_temp;
                fg_data_next = fg_data_temp;
                data_v_next = data_v_temp;
                bg_data_v_next = bg_data_v_temp;
                fg_data_v_next = fg_data_v_temp;
            else
                opTime2_next = [];
            end
        end

        %bg_data, 
        %fg_data, 
        %data_v, 
        %bg_data_v, 
        %fg_data_v
    
    end

    plotGraphs = false;
    
    if plotGraphs
        %% Plot smoothed
        figure;
        for k=1:nSizes+2

            if (k <= nSizes)

                subplot(nSizes+2,3,3*(k-1)+1);
                plot(opTime2,data_density(k,:)/(1),'Color',tColor(k,:));
                hold on;
                plot(opTime2,bg_density(k,:)/(1),'Color','black','LineStyle',':', 'LineWidth',2);
                title(['Diameter: ', num2str(diameters(k)), '\mum']);
                ylabel('#/m^3');
                xlabel('time')
                ylim_curr = ylim;
                ylim_curr = 1.0*[-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                subplot(nSizes+2,3,3*(k-1)+2);
                plot(opTime2,bg_density(k,:)/(1),'Color',tColor(k,:),'LineStyle',':');

                title(['Diameter: ', num2str(diameters(k)), '\mum']);
                ylabel('#/m^3');
                xlabel('time')
                ylim([ylim_curr]);

                subplot(nSizes+2,3,3*(k-1)+3);
                plot(opTime2,fg_density(k,:)/(1),'Color',tColor(k,:));

                title(['Diameter: ', num2str(diameters(k)), '\mum']);
                ylabel('#/m^3');
                xlabel('time')
                ylim([ylim_curr]);

            elseif k== nSizes+1
                subplot(nSizes+2,3,3*(k-1)+1);
                currentValid = ~all(isnan(data_temp),1);
                plot(opTime2,nansum(data_temp,1),'k');
                hold on;
                plot(opTime2,nansum(bg_data_temp,1),'Color','black','LineStyle',':', 'LineWidth',2);

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim_curr = ylim;
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                subplot(nSizes+2,3,3*(k-1)+2);
                plot(opTime2,nansum(bg_data_temp,1),'k');
                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);


                subplot(nSizes+2,3,3*(k-1)+3);
                plot(opTime2,nansum(fg_data_temp,1),'k');

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);
            elseif k == nSizes+2
                subplot(nSizes+2,3,3*(k-1)+1);
                currentValid = ~all(isnan(data_temp),1);
                plot(opTime2,nansum(data_v_temp,1),'k');
                hold on;
                plot(opTime2,nansum(bg_data_v_temp,1),'Color','black','LineStyle',':', 'LineWidth',2);

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim_curr = ylim;
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                subplot(nSizes+2,3,3*(k-1)+2);
                plot(opTime2,nansum(bg_data_v_temp,1),'k');
                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim(ylim_curr);


                subplot(nSizes+2,3,3*(k-1)+3);
                plot(opTime2,nansum(fg_data_v_temp,1),'k');

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim(ylim_curr);
            end
        end

        %% Plot annotations
        for k=1:nSizes+2
            for l=1:3
                subplot(nSizes+2,3,3*(k-1)+l);

                for m=1:size(eventNames,1)
                    if k==1
                        xline(eventTimes(m),'k--',eventNames{m,1});
                    else
                        xline(eventTimes(m),'k--');
                    end
                end
            end
        end
    end
    
end