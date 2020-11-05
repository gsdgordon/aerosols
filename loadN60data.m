function [times2, countssum, diameters2, shapes2] = loadN60data(folder, filePrefix, startTime, endTime)


    times = [];
    shapes =[];
    diameters = [];
    for k=1:20
       currentFile_vsp = [filePrefix, num2str(k), '.vsp'];
       filePath_vsp = fullfile(folder, currentFile_vsp);

       currentFile_txt = [filePrefix, num2str(k), '.txt'];
       filePath_txt = fullfile(folder, currentFile_txt);

       if (~exist(filePath_vsp))
           break;
       end

       if (k==1)
           settings_vsp = detectImportOptions('C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\20202810_01_N60\--3-Run1.vsp', 'FileType', 'delimitedtext'); %FIX can't just rely on this
           settings_txt = detectImportOptions(filePath_txt, 'FileType', 'text');
           settings_txt.DataLines = [1,17];
           settings_txt.VariableTypes{2} = 'char';
       end

       T_vsp = readtable(filePath_vsp, settings_vsp);
       T_txt = readtable(filePath_txt, settings_txt);

       dateStr = T_txt{5,1};
       dateStr = dateStr{1};
       timeStr = T_txt{5,2};
       timeStr = timeStr{1};

       frameLength = 60; %N60 frame
       cameraFrameRate = 30; %FPS
       frameEndTime = datetime([dateStr, ' ', timeStr]);
       frameStartTime = frameEndTime - duration(0,0,frameLength);

       if (size(T_vsp,1) > 6)
            %Data found
            subT = T_vsp(7:end,:);

            nVals = size(subT,1);

            for p=1:nVals
                currentRaw = subT{p,1};
                currentRaw2 = currentRaw{1};
                vals = strsplit(currentRaw2);

                currentFrame = str2double(vals{1});
                currentTime = duration(0,0,currentFrame/cameraFrameRate) + frameStartTime;
                currentPIdx = str2double(vals{2});
                currentArea = str2double(vals{3});
                currentDiameter = str2double(vals{4});
                currentShapeFactor = str2double(vals{5});

                times = [times;currentTime];
                diameters = [diameters; currentDiameter];
                shapes = [shapes; currentShapeFactor];
            end
       else
            %times = [times;frameStartTime; frameEndTime];
            %diameters = [diameters; 0; 0;];
            %shapes = [shapes; 0; 0;]; 
       end

    end

    times2 = startTime:seconds(1):endTime;
    times2 = times2';

    diameters2 = zeros(size(times2));
    shapes2 = zeros(size(times2));
    counts = zeros(size(times2));

    for k=1:size(times2,1)-1
        temp = times >= times2(k) & times < times2(k+1);

        for t = 1:size(temp,1)
            if (temp(t))
                times2 = [times2; times(t);];
                diameters2 = [diameters2; diameters(t);];
                shapes2 = [shapes2; shapes(t);];
                counts = [counts; 1];
            end
        end
    end

    [times2, I] = sort(times2, 'ascend');
    diameters2 = diameters2(I);
    shapes2 = shapes2(I);
    counts = counts(I);
    countssum = zeros(size(counts));
    countAvTime = seconds(30);

    startT = times2(1);
    startIdx = 1;
    for tIdx = 1:size(times2,1)
        currentT = times2(tIdx);

        if (currentT - startT) >= countAvTime
            cumsum = sum(counts(startIdx:tIdx));
            countssum(startIdx:tIdx) = cumsum;
            startT = currentT;
            startIdx = tIdx;
        end
    end


    figure;
    subplot(3,1,1);
    plot(times2, countssum);

    subplot(3,1,2);
    plot(times2, diameters2);

    subplot(3,1,3);
    plot(times2, shapes2);
end