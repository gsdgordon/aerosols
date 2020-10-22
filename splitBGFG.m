function [bg, fg] = splitBGFG(density, avSampleTime, validTimes)

    nSizes = size(density,1);
    bg = nan(size(density));
    
    for k=1:nSizes

        currentBG = density(k,validTimes);
        
        if all(isnan(currentBG))
            continue;
        end

        padSize = 100;
        currentBG = padarray(currentBG,[0,padSize],'replicate'); %FIX use median of last 3 values?
        currentBG = medfilt1(currentBG,15);

        sampleFreq = 1/avSampleTime; % in Hz
        cutoffFreq = 0.01; %in Hz


        currentBG = lowpass(currentBG,cutoffFreq,sampleFreq, 'ImpulseResponse','fir');


        currentBG = currentBG(padSize+1:end-padSize);
        bg(k,validTimes) = currentBG;
    end

    fg = density - bg;
    %fg_density(fg_density < 0) = 0; % Can't have negative signal
    
end