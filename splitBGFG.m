function [bg, fg] = splitBGFG(density, avSampleTime, validTimes)

    nSizes = size(density,1);
    bg = nan(size(density));
    
    for k=1:nSizes

        currentBG = density(k,validTimes);
        
        if all(isnan(currentBG))
            continue;
        end

        padSize = 100;
        
        firstVal = median(currentBG(1:5), 'omitnan');
        lastVal = median(currentBG(end-5:end), 'omitnan');
        
        medFiltSize = 15;
        currentBG = [firstVal, currentBG, lastVal];
        currentBG = padarray(currentBG,[0,padSize-1],'replicate'); %FIX use median of last 3 values?
        currentBG(isnan(currentBG)) = 0; % FIX should check beofre and after filter
        currentBG = medfilt1(currentBG,medFiltSize);

        sampleFreq = 1/avSampleTime; % in Hz
        cutoffFreq = 0.01; %in Hz

        currentBG = lowpass(currentBG,cutoffFreq,sampleFreq, 'ImpulseResponse','fir');
        currentBG = circshift(currentBG, -floor(medFiltSize/2)); % So background is taken as being before event.


        currentBG = currentBG(padSize+1:end-padSize);
        bg(k,validTimes) = currentBG;
    end

    fg = density - bg;
    %fg_density(fg_density < 0) = 0; % Can't have negative signal
    
end