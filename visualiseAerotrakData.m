% File showing example of how to process/plot data from AeroTrak, fit a
% lognormal distribution and creat a video
% Author: George Gordon
% Date: 30/09/2020


% Clear any old data
clc;
clear variables;
close all;

% Load data
Y = 2020;
M = 10;
D = 28;

P = 1;
P_N60 = 1;

datestring = sprintf('%0.4d%0.2d%0.2d',Y,M,D);

folder = ['C:\Users\george\OneDrive - The University of Nottingham\SAVE\',datestring,'\'];
%folder = ['/home/george/Desktop/'];

T = readtable(fullfile(folder, [datestring,'_aerotrak.xlsx']));
%T = readtable('/home/george/Desktop/20201030_aerotrak.xlsx');

opTime = T.DateAndTime;
location = T.Location;
sampleTime = T.SampleTime;
avSampleTime = mode(sampleTime); % assumes sample time is not changed during operation, but excludes partial samples

loadAnnotations = true;
if loadAnnotations
    if isnumeric(P)
        annotationFile = fullfile(folder,[datestring, '_patient', num2str(P), '.csv']);
    else
        annotationFile = fullfile(folder,[datestring, '_patient', P, '.csv']);
    end
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
    
    eventNames = T_annotation(validEvents,1);
    eventTimes_endo = T_annotation(validEvents,2);
    eventTimes_obscam = T_annotation(validEvents,3);
    eventTimes_aerotrak = T_annotation(validEvents,4);
    
    for k=1:size(eventNames,1)
        currentTime_dur =(eventTimes_aerotrak{k,1}); % FIX need to sort by date
        eventTimes(k,1) = datetime(Y,M,D,0,0,0) + currentTime_dur - seconds(avSampleTime); % Shift events, not aerotrak
    end
    
    % FIX need to sort by date
    startTime = eventTimes(1,1) - duration(0,5,0);
    endTime = eventTimes(end,1) + duration(0,5,0);
    
    aeroOffsetTime = 0;
else
    % Sync times from different clocks as per the video
    obsCamTime_endo = datetime(Y,M,D,11,27,35);
    endoscopeTime = datetime(Y,M,D,11,41,30);

    obsCamTime_aerotrak = datetime(Y,M,D,11,27,39);
    aerotrakTime = datetime(Y,M,D,11,39,04) + seconds(avSampleTime); %Aerotrak time is the time at the start of the sample

    aeroOffsetTime = aerotrakTime - obsCamTime_aerotrak;
    endoOffsetTime = endoscopeTime - obsCamTime_endo;

    opTime = opTime - aeroOffsetTime;

    startTime = datetime(Y,M,D,11,35,00);
    endTime = datetime(Y,M,D,11,50,00);
end


tValid = isbetween(opTime,startTime,endTime);
tValid = tValid & strcmpi(location,'Location01');

T = T(tValid,:);
T = T(1:end-1,:); % remove last count as it is likely partial

opTime2 = T.DateAndTime - aeroOffsetTime;
airVol = T.Volume_L_;

loadN60data_var = true;
if loadN60data_var
    fileListN60_par = dir(folder);
    filterFun = @(x) regexpi(x, [datestring, '_', sprintf('%02d',P_N60), '_N60.*']);
    temp = cellfun(filterFun, {fileListN60_par.name}, 'UniformOutput', false); 
    fileListN60_par = fileListN60_par(~cellfun(@isempty,temp));
        
    folderN60 = [folder, fileListN60_par(1).name];
    
    fileListN60 = dir(folderN60);
    
    filterFun = @(x) regexpi(x, '--([0-9])+.*', 'match', 'tokens');
    [temp, temp2] = cellfun(filterFun, {fileListN60.name}, 'UniformOutput', false); 
    temp2 = temp2';
    
    for k = 1:size(temp2,1)
        currentNum = temp2(k);
        currentNum = currentNum{1};
            
        if ~isempty(currentNum)
            currentNum = currentNum{1};
            currentNum = currentNum{1};

            prefixList(k) = str2num(currentNum);
        end
    end
    
    prefixNum = mode(prefixList);

    prefix = ['--', sprintf('%d',prefixNum), '-Run'];

    AT_diff = eventTimes_aerotrak{1,1} - eventTimes_obscam{1,1};

    obscamTime_N60 = datetime(2020,10,28,9,37,35);
    N60_time = datetime(2020,10,28,8,48,56);
    N60_diff = N60_time - obscamTime_N60;

    N60_diff_AT = N60_diff - AT_diff;

    % INput times
    if (~isempty(opTime2))
        startN60 = opTime2(1)+N60_diff_AT;
        endN60  = opTime2(end)+N60_diff_AT;
    else
        startN60 = startTime - aeroOffsetTime + N60_diff_AT;
        endN60  = endTime - aeroOffsetTime + N60_diff_AT;
    end
    
    useEvent = true;
    if useEvent
        eventToUse = 'Cough';
        idx = find(strcmpi(table2cell(eventNames),eventToUse));
        idx = idx(1);
        startN60 = eventTimes(idx) + N60_diff_AT;
        endN60 = startN60 + seconds(60);
    end
        
    [N60timesall, N60counts, N60diameters, N60shapes] = loadN60data(folderN60, prefix, startN60, endN60);
    N60timesall = N60timesall - N60_diff_AT;
    
    N60_area = 12606*7181; %0.58 magnification, in micronrs
    N60_depth = 50000;
    N60_vol = N60_area/(1e6^2) * N60_depth/(1e6);
    
    if useEvent
        save([folder, '/N60_agg_', num2str(P_N60), '_', eventToUse, '.mat'], 'N60diameters', 'N60shapes');
    else
        save([folder, '/N60_agg_', num2str(P_N60), '.mat'], 'N60diameters', 'N60shapes');
    end
    
    %"Max. Probe volume (cu mm)"	4526.4189
end

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


%% Plot smoothed
figure;
for k=1:nSizes+2

    if (k <= nSizes)

        subplot(nSizes+2,3,3*(k-1)+1);
        plot(opTime2,data_density(k,:)/(1-correctionVal(k)),'Color',tColor(k,:));
        hold on;
        plot(opTime2,bg_density(k,:)/(1-correctionVal(k)),'Color','black','LineStyle',':', 'LineWidth',2);
        title(['Diameter: ', num2str(diameters(k)), '\mum']);
        ylabel('#/m^3');
        xlabel('time')
        ylim_curr = ylim;
        ylim_curr = 1.0*[-1*max(ylim_curr),max(ylim_curr)];
        ylim(ylim_curr);

        subplot(nSizes+2,3,3*(k-1)+2);
        plot(opTime2,bg_density(k,:)/(1-correctionVal(k)),'Color',tColor(k,:),'LineStyle',':');

        title(['Diameter: ', num2str(diameters(k)), '\mum']);
        ylabel('#/m^3');
        xlabel('time')
        ylim([ylim_curr]);

        subplot(nSizes+2,3,3*(k-1)+3);
        plot(opTime2,fg_density(k,:)/(1-correctionVal(k)),'Color',tColor(k,:));

        title(['Diameter: ', num2str(diameters(k)), '\mum']);
        ylabel('#/m^3');
        xlabel('time')
        ylim([ylim_curr]);
        
    elseif k== nSizes+1
        subplot(nSizes+2,3,3*(k-1)+1);
        currentValid = ~all(isnan(data),1);
        plot(opTime2,nansum(data./repmat(1-correctionVal',1,size(data,2)),1),'k');
        hold on;
        plot(opTime2,nansum(bg_data,1),'Color','black','LineStyle',':', 'LineWidth',2);
        
        title(['Total #']);
        ylabel('#/m^3');
        xlabel('time');
        ylim_curr = ylim;
        ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
        ylim(ylim_curr);

        subplot(nSizes+2,3,3*(k-1)+2);
        plot(opTime2,nansum(bg_data./repmat(1-correctionVal',1,size(data,2)),1),'k');
        title(['Total #']);
        ylabel('#/m^3');
        xlabel('time');
        ylim(ylim_curr);


        subplot(nSizes+2,3,3*(k-1)+3);
        plot(opTime2,nansum(fg_data./repmat(1-correctionVal',1,size(data,2)),1),'k');

        title(['Total #']);
        ylabel('#/m^3');
        xlabel('time');
        ylim(ylim_curr);
    elseif k == nSizes+2
        subplot(nSizes+2,3,3*(k-1)+1);
        currentValid = ~all(isnan(data),1);
        plot(opTime2,nansum(data_v./repmat(1-correctionVal',1,size(data,2)),1),'k');
        hold on;
        plot(opTime2,nansum(bg_data_v./repmat(1-correctionVal',1,size(data,2)),1),'Color','black','LineStyle',':', 'LineWidth',2);
        
        title(['Total vol']);
        ylabel('vol/m^3');
        xlabel('time');
        ylim_curr = ylim;
        ylim_curr = [-1*max(ylim_curr),max(ylim_curr)];
        ylim(ylim_curr);

        subplot(nSizes+2,3,3*(k-1)+2);
        plot(opTime2,nansum(bg_data_v./repmat(1-correctionVal',1,size(data,2)),1),'k');
        title(['Total vol']);
        ylabel('vol/m^3');
        xlabel('time');
        ylim(ylim_curr);


        subplot(nSizes+2,3,3*(k-1)+3);
        plot(opTime2,nansum(fg_data_v./repmat(1-correctionVal',1,size(data,2)),1),'k');

        title(['Total vol']);
        ylabel('vol/m^3');
        xlabel('time');
        ylim(ylim_curr);
    end
end

for k=1:nSizes+2
    for l=1:3
        subplot(nSizes+2,3,3*(k-1)+l);
        % Patient 1
%         xline(datetime(2020,10,30,10,11,22),'k--','e');
%         xline(datetime(2020,10,30,10,11,29),'k--','m');
%         xline(datetime(2020,10,30,9,59,39),'k--','sp');
%         xline(datetime(2020,10,30,10,00,46),'k--','mv');
%         xline(datetime(2020,10,30,10,03,39),'k--','i');

        
        %Patient 2
%         xline(datetime(2020,10,30,12,00,44),'k--','e');
%         xline(datetime(2020,10,30,12,00,50),'k--','m');
%         xline(datetime(2020,10,30,11,48,00),'k--','sp');
%         xline(datetime(2020,10,30,11,50,35),'k--','mv');
%         xline(datetime(2020,10,30,11,52,15),'k--','i');
        
        % Patient 3
%         xline(datetime(2020,10,30,13,00,17),'k--','e');
%         xline(datetime(2020,10,30,12,45,28),'k--','sp');
%         xline(datetime(2020,10,30,12,47,22),'k--','mv');
%         xline(datetime(2020,10,30,13,00,38),'k--','m');
%         xline(datetime(2020,10,30,13,01,00),'k--','mv');
%         xline(datetime(2020,10,30,12,49,46),'k--','i');

        % Patient 4 (no mask)
%         xline(datetime(2020,10,22,16,10,28),'k--','e');
%         xline(datetime(2020,10,22,15,55,03),'k--','sp');
%         xline(datetime(2020,10,22,15,56,05),'k--','mv');
%         xline(datetime(2020,10,22,16,01,10),'k--','i');
        
        % Patient 5 (no mask)
%         xline(datetime(2020,10,22,17,10,59),'k--','e');
%         xline(datetime(2020,10,22,16,54,24),'k--','sp');
%         xline(datetime(2020,10,22,16,56,06),'k--','mv');
%         xline(datetime(2020,10,22,17,03,37),'k--','i');
        
        % Patient 6 (no mask)
%         xline(datetime(2020,10,22,18,05,42),'k--','e');
%         xline(datetime(2020,10,22,17,50,25),'k--','sp');
%         xline(datetime(2020,10,22,17,52,20),'k--','mv');
%         xline(datetime(2020,10,22,17,54,55),'k--','i');

        for m=1:size(eventNames,1)
            if k==1
                xline(eventTimes(m),'k--',eventNames{m,1});
            else
                xline(eventTimes(m),'k--');
            end
        end
    end
end

if loadN60data_var
    figure;
    %N60timesall(567) = N60timesall(567)+seconds(0.01);
    N60_countsResamp = interp1(N60timesall, N60counts, opTime2);
    largeParts = data_density(6,:)/(1-correctionVal(6));
    largeParts = largeParts.';
    
    [corrCurve, lags] = xcorr(largeParts, N60_countsResamp);
    
    lags = lags.' .* avSampleTime;
    winWidth = 180;
    valCorr = lags >= -1*winWidth & lags <= winWidth;
    
    corrCurve = corrCurve(valCorr);
    lags = lags(valCorr);
    
    temp = corrcoef(largeParts, N60_countsResamp);
    normVal = temp(1,2);
    
    corrCurve = corrCurve/corrCurve(lags == 0) * normVal;
    
    plot(lags, corrCurve);
    xlabel('offset in seconds');
    ylabel('correlation');
    
    save([folder, '/N60_corr_', num2str(P_N60), '.mat'], 'lags', 'corrCurve');
    
    figure;
    subplot(7,1,1);
    plot(N60timesall, N60counts);
    title('N60 counts');
    xlabel('time');
    ylabel('count');

    subplot(7,1,2);
    plot(N60timesall, N60diameters, 'r');
    title('N60 diameters');
    xlabel('time');
    ylabel('diameter (\mum)');

    subplot(7,1,3);
    plot(N60timesall, N60shapes,'b');
    title('N60 shapes');
    xlabel('time');
    ylabel('sphericity');
    
    for k=5:nSizes+2

        if (k <= nSizes)

            subplot(7,1,k-1);
            plot(opTime2,data_density(k,:)/(1-correctionVal(k)),'Color',tColor(k,:));
            title(['AeroTrak Diameter: ', num2str(diameters(k)), '\mum']);
            ylabel('#/m^3');
            xlabel('time')


        elseif k== nSizes+1
            subplot(7,1,k-1);
            currentValid = ~all(isnan(data),1);
            plot(opTime2,nansum(data./repmat(1-correctionVal',1,size(data,2)),1),'k');


            title(['AeroTrak Total #']);
            ylabel('#/m^3');
            xlabel('time');
            

        elseif k == nSizes+2
            subplot(7,1,k-1)
            currentValid = ~all(isnan(data),1);
            plot(opTime2,nansum(data_v./repmat(1-correctionVal',1,size(data,2)),1),'k');

            title(['AeroTrak Total vol']);
            ylabel('vol/m^3');
            xlabel('time');
        end
    end
    
    for k=1:nSizes+2
        if k <= 3
            subplot(7,1,k);
        else
        	subplot(7,1,k-1);
        end
        
        for m=1:size(eventNames,1)
            if k==1
                xline(eventTimes(m),'k--',eventNames{m,1});
            else
                xline(eventTimes(m),'k--');
            end
        end
    end
    
end

%save('cond4.mat','data', 'diameters', 'opTime2','avSampleTime');

% figure;
% %Try to fit lognormal distribution to BG
% initPop = [];
% for k=1:size(bg_density,2)
%     
%     densities = bg_density(:,k);
%     
%     if any(isnan(densities))
%         A_bg(k) = NaN;
%         mu_bg(k) = NaN;
%         sigma_bg(k) = NaN;
%     else
%         if k>1
%             initPop = [mu_bg(1:k-1)',sigma_bg(1:k-1)']; % Avoid nans!
%             validRows = ~isnan(initPop(:,1)) & ~isnan(initPop(:,2));
%             initPop = initPop(validRows,:);
%         end
%         
%         [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters, densities','fitType','counts','initPop',initPop, 'mu_LB', log(0.01), 'mu_UB', log(1), 'sig_LB', 0.01, 'sig_UB', 1.5, 'mu_mu', log(0.2), 'mu_sig',5);
%         
%         A_bg(k) = A_t;
%         mu_bg(k) = mu_t;
%         sigma_bg(k) = sigma_t;
%         likelihood_bg(k) = l_t;
%         
%         if mod(k,10) == 0
%             disp(['Done ', num2str(k), '/', num2str(size(data,2))]);
%         end
% 
%     end
% 
% end
% 
% figure;
% subplot(4,1,1);
% plot(A_bg);
% title('BG Aerosol number');
% 
% subplot(4,1,2);
% plot(mu_bg);
% title('BG mu')
% 
% subplot(4,1,3);
% plot(sigma_bg);
% title('BG sigma')
% 
% subplot(4,1,4);
% plot(likelihood_bg);
% title('BG log likelihood');
% 
% mean_mu = median(mu_bg); % Do unbiased estimate!
% mean_sig = median(sigma_bg);
% 
% initPop = [];

clipNegatives = true;

figure;
%Try to fit lognormal distribution to foreground
startK = 1;
for k=startK:size(data,2)
    
    currentLogBinSizes = log(diameters(2:end)) - log(diameters(1:end-1));
    currentLogBinSizes = currentLogBinSizes';
    currentData = fg_data(1:end,k);
    currentData_raw = currentData;
    currentDiameters = diameters';
    
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
    
    %%%%%
    
    if any(isnan(currentData))
        A(k) = NaN;
        mu(k) = NaN;
        sigma(k) = NaN;
    else
        if k>startK
            initPop = [mu(1:k-1)',sigma(1:k-1)']; % Avoid nans!
            validRows = ~isnan(initPop(:,1)) & ~isnan(initPop(:,2));
            initPop = initPop(validRows,:);
        else
            initPop = [];
        end
        
        [~, mu_t, sigma_t, l_t] = fitAerosolDist(currentDiameters.', (currentData./currentLogBinSizes)','fitType','counts','initPop',initPop, 'mu_LB', log(0.01), 'mu_UB', log(50), 'sig_LB', 0.1, 'sig_UB', 50);
        
        
        A(k) = A_t;
        mu(k) = mu_t;
        sigma(k) = sigma_t;
        likelihood_fg(k) = l_t;

    
        if mod(k,10) == 0
            disp(['Done ', num2str(k), '/', num2str(size(data,2))]);
        end

    end

end

isValid = A > 4e4;

subplot(4,1,1)
plot(opTime2, A.*isValid);
title('Relative aerosol volume');

subplot(4,1,2);
plot(opTime2, exp(mu).*isValid);
title('\mu (mean aerosol diameter in microns)');

subplot(4,1,3);
plot(opTime2, sigma.*isValid);
title('\sigma (standard devation of particle diameters in microns)');

subplot(4,1,4);
plot(opTime2, (likelihood_fg).*isValid);
title('log likelihood');




data_tot = sum(data_v,2);


tColor = lines(nSizes);

allFig = figure;
sumFig = figure;

figure(allFig);
for k=1:nSizes+1
    
    if (k <= nSizes)
        subplot(nSizes,1,k);
        plot(opTime2,data(:,k),'Color',tColor(k,:));
        title(['Diameter: ', num2str(diameters(k)), '\mum']);
        ylabel('#/m^3');
        xlabel('time')
       
        if k==1
            %yline(10200,'k:', 'ISO5')
            yline(102000,'k:', 'ISO6')
        end
        
        if k==2
            %yline(3520,'k:', 'ISO5')
            yline(35200,'k:', 'ISO6')
            yline(352000,'k:', 'ISO7')
        end
        
        if k==3
            %yline(832,'k:', 'ISO5')
            %yline(8320,'k:', 'ISO6')
            yline(83200,'k:', 'ISO7')
        end
        
        if k==5
            %yline(29,'k:', 'ISO5')
            %yline(293,'k:', 'ISO6')
            yline(2930,'k:', 'ISO7')
        end
        
        ylim_array(k,1:2) = ylim;
    else
        figure(sumFig)
        plot(opTime2,data_tot,'Color','k');
    end
    
    %xline(datetime(2020,8,7,14,33,00),'r--','attached to bed');


end


% Make a video suitable for overlaying on clinical videos
makeVideo = false;

if makeVideo
    fps = 30;
    window = 60; % in seconds

    d = seconds(endTime - startTime);
    nFrames = d*fps;

    v = VideoWriter('test', 'FrameRate', fps);
    open(v);
    
    for frameIdx = 1:nFrames
        
        currentStart = startTime + duration(0,0,frameIdx/fps - window);
        currentEnd = currentStart + duration(0,0,window);
    
        for k=1:nSizes

            subplot(nSizes,1,k);
            if (frameIdx == 1)
                
                plot(opTime2,data(:,k),'Color',tColor(k,:));
                %title(['Diameter: ', num2str(diameters(k)), '\mum']);
                ylabel('v/m^3');
                %xlabel("time")
                
                set(gcf,'Position',[100,100,200,600]);
            end
            axis tight;
            ylim(ylim_array(k,1:2));
            xlim([currentStart,currentEnd]);
            
            ax = gca;
            set(ax(1),'XTickLabel','')

        end
        
        frame = getframe(gcf);
        writeVideo(v,frame);
        

        
        if mod(frameIdx,fps*60) == 0
            disp(['Frame ', num2str(frameIdx)]);
            pause(1/fps);
        end
        
    end
    
    close(v);

    
end

