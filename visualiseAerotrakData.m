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
M = 9;
D = 30;

datestring = sprintf('%0.4d%0.2d%0.2d',Y,M,D);

T = readtable([datestring,'_aerotrak.xlsx']);

opTime = T.DateAndTime;
location = T.Location;
sampleTime = T.SampleTime;
avSampleTime = mode(sampleTime); % assumes sample time is not changed during operation, but excludes partial samples


% Sync times from different clocks as per the video
obsCamTime_endo = datetime(Y,M,D,11,09,55);
endoscopeTime = datetime(Y,M,D,11,24,20);

obsCamTime_aerotrak = datetime(Y,M,D,11,10,01);
aerotrakTime = datetime(Y,M,D,11,22,33) + seconds(avSampleTime); %Aerotrak time is the time at the start of the sample

aeroOffsetTime = aerotrakTime - obsCamTime_aerotrak;
endoOffsetTime = endoscopeTime - obsCamTime_endo;

opTime = opTime - aeroOffsetTime;

startTime = datetime(Y,M,D,11,20,00);
endTime = datetime(Y,M,D,12,50,00);

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
diameters_av = (diameters(1:6)+diameters(2:7))/2; % Mean diameter in each bin


% The actual counts in each bin
data = [T.Ch1Diff___,T.Ch2Diff___,T.Ch3Diff___,T.Ch4Diff___,T.Ch5Diff___,T.Ch6Diff___];
data = data./repmat(airVol,1,size(data,2)) * 1000; %Because volume is liters so times by 1000 to get to m^3


% Convert counts to volumes
vols = 4/3*pi*(diameters_av/2).^3 * (1e-6)^3;

bin_sizes = diameters(2:7) - diameters(1:6);
log_bin_sizes = log(diameters(2:7)) - log(diameters(1:6));

data_v = data .* repmat(vols,size(data,1),1);

% Densities so that a probability density approach can be used
data_v_density = data_v ./ repmat(log_bin_sizes,size(data,1),1); %Try using log binsizes
data_density = data ./ repmat(log_bin_sizes,size(data,1),1);
%data_v_density = data_v ./ repmat(bin_sizes,size(data,1),1); %Try using linear binsizes
%data_density = data ./ repmat(bin_sizes,size(data,1),1);

nSizes = size(data,2);
tColor = lines(nSizes);

% Get a background reading
%bgStartTime = datetime(2020,9,30,11,21,00);
%bgEndTime = datetime(2020,9,30,11,23,00);
bgStartTime = startTime;
bgEndTime = endTime;
bgValid = isbetween(opTime2,bgStartTime,bgEndTime);
bg_density_raw = data_density(bgValid,:);


for k=1:nSizes
  
    currentBG = bg_density_raw(:,k);
    
    padSize = 100;
    currentBG = padarray(currentBG,padSize,'replicate');
    currentBG = medfilt1(currentBG,15);
    
    sampleFreq = 1/avSampleTime; % in Hz
    cutoffFreq = 0.01; %in Hz
    

    currentBG = lowpass(currentBG,cutoffFreq,sampleFreq, 'ImpulseResponse','fir');
    
    
    currentBG = currentBG(padSize+1:end-padSize);
    bg_density_filt(:,k) = currentBG;
end

fg_density = bg_density_raw - bg_density_filt;
fg_density(fg_density < 0) = 0; % Can't have negative signal

% Plot smoothed
figure;
for k=1:nSizes
    
    subplot(nSizes,1,k);
    plot(opTime2,bg_density_raw(:,k),'Color',tColor(k,:));
    hold on;
    plot(opTime2,bg_density_filt(:,k),'Color','black','LineStyle',':', 'LineWidth',2);

    title(['Diameter: ', num2str(diameters(k)), '\mum']);
    ylabel('#/m^3');
    xlabel('time')

end

centreTime_idx = 76;
resampleHalfWidthTime = 180;
startTime_idx = ceil(centreTime_idx - resampleHalfWidthTime/avSampleTime);
endTime_idx = floor(centreTime_idx + resampleHalfWidthTime/avSampleTime);
data_samp = data(startTime_idx:endTime_idx,:);
data_t = ((startTime_idx:1:endTime_idx) - centreTime_idx)*avSampleTime;


new_t = -resampleHalfWidthTime:resampleHalfWidthTime;
new_data = zeros(size(new_t,2), size(data_samp,2));

count = 1;
for tIdx = 1:size(new_t,2)
    if count <= size(data_t,2)
        if data_t(count) == new_t(tIdx)
            new_data(tIdx,:) = data_samp(count,:);
        	count = count + 1;
        else
            new_data(tIdx,:) = NaN;
        end
    else
        new_data(tIdx,:) = NaN;
    end
end

new_data = new_data';

% Plot FG
figure;
for k=1:nSizes
    
    subplot(nSizes,1,k);
    plot(opTime2,fg_density(:,k),'Color',tColor(k,:));

    title(['Diameter: ', num2str(diameters(k)), '\mum']);
    ylabel('#/m^3');
    xlabel('time')

end

figure;
%Try to fit lognormal distribution to BG
initPop = [];
for k=1:size(bg_density_filt,1)
    
    densities = bg_density_filt(k,:);
    
    if any(isnan(densities))
        A_bg(k) = NaN;
        mu_bg(k) = NaN;
        sigma_bg(k) = NaN;
    else
        if k>1
            initPop = [mu_bg(1:k-1)',sigma_bg(1:k-1)']; % Avoid nans!
            validRows = ~isnan(initPop(:,1)) & ~isnan(initPop(:,2));
            initPop = initPop(validRows,:);
        end
        
        [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters, densities,'fitType','counts','initPop',initPop, 'mu_LB', log(0.01), 'mu_UB', log(1), 'sig_LB', 0.01, 'sig_UB', 1.5, 'mu_mu', log(0.2), 'mu_sig',5);
        
        A_bg(k) = A_t;
        mu_bg(k) = mu_t;
        sigma_bg(k) = sigma_t;
        likelihood_bg(k) = l_t;
        
        if mod(k,10) == 0
            disp(['Done ', num2str(k), '/', num2str(size(data,1))]);
        end

    end

end

figure;
subplot(4,1,1);
plot(A_bg);
title('BG Aerosol number');

subplot(4,1,2);
plot(mu_bg);
title('BG mu')

subplot(4,1,3);
plot(sigma_bg);
title('BG sigma')

subplot(4,1,4);
plot(likelihood_bg);
title('BG log likelihood');

mean_mu = median(mu_bg); % Do unbiased estimate!
mean_sig = median(sigma_bg);

initPop = [];

figure;
%Try to fit lognormal distribution
startK = 285;
for k=startK:size(data,1)
    
    densities = fg_density(k,:);
    %densities_v = data_v_density(k,:);
    
    if any(isnan(densities))
        A(k) = NaN;
        mu(k) = NaN;
        sigma(k) = NaN;
    else
        if k>startK
            initPop = [mu(1:k-1)',sigma(1:k-1)']; % Avoid nans!
            validRows = ~isnan(initPop(:,1)) & ~isnan(initPop(:,2));
            initPop = initPop(validRows,:);
        end
        
        [A_t, mu_t, sigma_t, l_t] = fitAerosolDist(diameters, densities,'fitType','counts','initPop',initPop, 'mu_LB', log(0.01), 'mu_UB', log(50), 'sig_LB', 0.1, 'sig_UB', 50);
        
        
        A(k) = A_t;
        mu(k) = mu_t;
        sigma(k) = sigma_t;
        likelihood_fg(k) = l_t;

    
        if mod(k,10) == 0
            disp(['Done ', num2str(k), '/', num2str(size(data,1))]);
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

