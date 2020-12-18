function [data, opTime2, eventTimes, eventNames, avSampleTime, bg_data, fg_data, data_v, bg_data_v, fg_data_v, diameters, diameters_av] = loadAnnotatedData(Y,M,D,P)

    datestring = sprintf('%0.4d%0.2d%0.2d',Y,M,D);

    folder = ['C:\Users\george\OneDrive - The University of Nottingham\SAVE\',datestring,'\'];

    T = readtable(fullfile(folder, [datestring,'_aerotrak.xlsx']));
    %T = readtable('/home/george/Desktop/20201030_aerotrak.xlsx');

    opTime = T.DateAndTime;
    location = T.Location;
    sampleTime = T.SampleTime;
    avSampleTime = mode(sampleTime); % assumes sample time is not changed during operation, but excludes partial samples

    annotationFile = fullfile(folder,[datestring, '_patient', num2str(P), '.csv']);
    T_annotation = readtable(annotationFile, 'ReadVariableNames', false, 'HeaderLines', 0); % FIX ignores first row

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
        return;
    end

    eventNames = T_annotation(validEvents,1);
    %eventTimes_endo = T_annotation(validEvents,2);
    %eventTimes_obscam = T_annotation(validEvents,3);
    eventTimes_aerotrak = T_annotation(validEvents,4);

    for k=1:size(eventNames,1)
        currentTime_dur =(eventTimes_aerotrak{k,1}); % FIX need to sort by date
        eventTimes(k,1) = datetime(Y,M,D,0,0,0) + currentTime_dur - seconds(avSampleTime); % Shift events, not aerotrak
    end
    
    [eventTimes, sortOrder] = sort(eventTimes,'ascend');
    eventNames = eventNames(sortOrder,1);

    % FIX need to sort by date
    startTime = eventTimes(1,1) - duration(0,5,0);
    endTime = eventTimes(end,1) + duration(0,5,0);

    aeroOffsetTime = 0;

    tValid = isbetween(opTime,startTime,endTime);
    tValid = tValid & strcmpi(location,'Location01');

    T = T(tValid,:);
    T = T(1:end-1,:); % remove last count as it is likely partial

    opTime2 = T.DateAndTime - aeroOffsetTime;
    airVol = T.Volume_L_;

    % Now get the particle 'diameters' representing the edges of the counting
    % bins
    maxDiameter = 25; % in microns - this is from the spec sheet
    diameters = [T.Ch1Size__m_(1),T.Ch2Size__m_(1),T.Ch3Size__m_(1),T.Ch4Size__m_(1),T.Ch5Size__m_(1),T.Ch6Size__m_(1), maxDiameter];
    diameters_av = (diameters(1:end-1)+diameters(2:end))/2; % Mean diameter in each bin

    % The actual counts in each bin
    data = [T.Ch1Diff___,T.Ch2Diff___,T.Ch3Diff___,T.Ch4Diff___,T.Ch5Diff___,T.Ch6Diff___];
    data = data./repmat(airVol,1,size(data,2)) * 1000; %Because volume is liters so times by 1000 to get to m^3
    data = data';

    [bg_data, fg_data] = splitBGFG(data, avSampleTime, true(1,size(data,2)));

    % Convert counts to volumes
    vols = 4/3*pi*(diameters_av/2).^3 * (1e-6)^3;

    bin_sizes = diameters(2:7) - diameters(1:6);
    log_bin_sizes = log(diameters(2:7)) - log(diameters(1:6));

    data_v = data .* repmat(vols',1,size(data,2));
    bg_data_v = bg_data .* repmat(vols',1,size(data,2));
    fg_data_v = fg_data .* repmat(vols',1,size(data,2));

    % Densities so that a probability density approach can be used
    data_v_density = data_v ./ repmat(log_bin_sizes',1,size(data,2)); %Try using log binsizes
    data_density = data ./ repmat(log_bin_sizes',1,size(data,2));
    %data_v_density = data_v ./ repmat(bin_sizes,size(data,1),1); %Try using linear binsizes
    %data_density = data ./ repmat(bin_sizes,size(data,1),1);

    bg_density = bg_data ./ repmat(log_bin_sizes',1,size(data,2));
    fg_density = fg_data ./ repmat(log_bin_sizes',1,size(data,2));

    nSizes = size(data,1);
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
    
    data = data ./ repmat(1-correctionVal',1,size(data,2));
    bg_data = bg_data ./ repmat(1-correctionVal',1,size(data,2));
    fg_data = fg_data ./ repmat(1-correctionVal',1,size(data,2));
    data_v = data_v ./ repmat(1-correctionVal',1,size(data,2));
    bg_data_v = bg_data_v ./ repmat(1-correctionVal',1,size(data,2));
    fg_data_v = fg_data_v ./ repmat(1-correctionVal',1,size(data,2));
   
    %bg_data, 
    %fg_data, 
    %data_v, 
    %bg_data_v, 
    %fg_data_v

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
                currentValid = ~all(isnan(data),1);
                plot(opTime2,nansum(data,1),'k');
                hold on;
                plot(opTime2,nansum(bg_data,1),'Color','black','LineStyle',':', 'LineWidth',2);

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim_curr = ylim;
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                subplot(nSizes+2,3,3*(k-1)+2);
                plot(opTime2,nansum(bg_data,1),'k');
                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);


                subplot(nSizes+2,3,3*(k-1)+3);
                plot(opTime2,nansum(fg_data,1),'k');

                title(['Total #']);
                ylabel('#/m^3');
                xlabel('time');
                ylim(ylim_curr);
            elseif k == nSizes+2
                subplot(nSizes+2,3,3*(k-1)+1);
                currentValid = ~all(isnan(data),1);
                plot(opTime2,nansum(data_v,1),'k');
                hold on;
                plot(opTime2,nansum(bg_data_v,1),'Color','black','LineStyle',':', 'LineWidth',2);

                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim_curr = ylim;
                ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
                ylim(ylim_curr);

                subplot(nSizes+2,3,3*(k-1)+2);
                plot(opTime2,nansum(bg_data_v,1),'k');
                title(['Total vol']);
                ylabel('vol/m^3');
                xlabel('time');
                ylim(ylim_curr);


                subplot(nSizes+2,3,3*(k-1)+3);
                plot(opTime2,nansum(fg_data_v,1),'k');

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