clc;
clear variables;
close all;

folder = 'C:\Users\george\OneDrive - The University of Nottingham\SAVE\01-27-2021_results\Summary_LT5\';
load([folder, 'resultsTable_lt5.mat']);

% Plot var importances
for k = 1:8
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    % Plot table
    if k == 1
        varName = 'predictorImportance';
        colScale = [0,15e8];
    elseif k == 2
        varName = 'predictorImportance_v';
        colScale = [0, 1e-22];
    elseif k == 3
        varName = 'predictorImportance_mu';
        colScale = [0, 0.7];
    elseif k == 4
        varName = 'predictorImportance_sig';
        colScale = [0,27];
    elseif k == 5
        varName = 'predictorImportance_perm';
        colScale = [-1,1];
    elseif k == 6
        varName = 'predictorImportance_perm_v';
        colScale = [-1, 1];
    elseif k == 7
        varName = 'predictorImportance_perm_mu';
        colScale = [-1, 1];
    elseif k == 8
        varName = 'predictorImportance_perm_sig';
        colScale = [-1,1];
    end
    eval(['temp = resultsTable.' varName ';']);
    temp = temp(3:end);
    temp2 = cell2mat(temp);

    %temp2(temp2 < 0)  = 0;
    
    % Lower
    lgi = temp2(1:10,:);
    lgi_eventvalid = logical([1;1;1;0;0;0;0;1;1;01;]);
    lgi = lgi(lgi_eventvalid,:);
    lgi_varvalid = logical([1;1;1;1;1;1;1;0;0;0;1;1;1;0;0]);
    lgi = lgi(:,lgi_varvalid);
    
    subplot(2,2,1);
    imagesc((lgi), colScale);
    %imagesc((temp2));
    h = gca;
    %axis equal;

    varNames = {'age', 'sex', 'bmi', 'sedation', 'analTone', 'co2water', 'isSmoker', 'patientMask', 'roomType', 'ugiRoute', 'divertDisease', 'looping', 'discomfort', 'hiatusHernia', 'suctioning'};

    varNames_temp = varNames(lgi_varvalid.');
    h.XTickLabel = varNames_temp;
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = 'none';
    xticks(1:size(varNames_temp,2));

    eventNames = resultsTable.label;
    eventNames = eventNames(3:end);
    eventNames = eventNames(1:10);
    eventNames_temp = eventNames(lgi_eventvalid);
    h.YTickLabel = eventNames_temp;
    yticks(1:size(eventNames_temp,1));
    colorbar;
    
    title(varName);
    
    subplot(2,2,2);
    maxImp = max(lgi,[],1);
    meanImp = mean(lgi,1);
    
    bar([maxImp; meanImp]');
    h = gca;
    h.XTickLabel = varNames_temp;
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = 'none';
    xticks(1:size(varNames_temp,2));
    ylim(colScale)
    
    legend('Max', 'Mean');
    
    % Upper
    ugi = temp2(11:end,:);
    ugi_eventvalid = logical([0;1;1;1;1;1;1;0;0;1;1]);
    ugi = ugi(ugi_eventvalid,:);
    ugi_varvalid = logical([1;1;1;1;0;0;1;0;0;0;0;0;1;1;0]);
    ugi = ugi(:,ugi_varvalid);
    
    subplot(2,2,3);
    imagesc((ugi), colScale);
    %imagesc((temp2));
    h = gca;
    %axis equal;

    varNames = {'age', 'sex', 'bmi', 'sedation', 'analTone', 'co2water', 'isSmoker', 'patientMask', 'roomType', 'ugiRoute', 'divertDisease', 'looping', 'discomfort', 'hiatusHernia', 'suctioning'};

    varNames_temp = varNames(ugi_varvalid.');
    h.XTickLabel = varNames_temp;
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = 'none';
    xticks(1:size(varNames_temp,2));

    eventNames = resultsTable.label;
    eventNames = eventNames(3:end);
    eventNames = eventNames(11:end);
    eventNames_temp = eventNames(ugi_eventvalid);
    h.YTickLabel = eventNames_temp;
    yticks(1:size(eventNames_temp,1));
    colorbar;
    
    title(varName);
    
    subplot(2,2,4);
    maxImp = max(ugi,[],1);
    meanImp = mean(ugi,1);
    
    bar([maxImp; meanImp]');
    h = gca;
    h.XTickLabel = varNames_temp;
    h.XTickLabelRotation = 45;
    h.TickLabelInterpreter = 'none';
    xticks(1:size(varNames_temp,2));
    ylim(colScale);
     
    legend('Max', 'Mean');
    
    saveFigName = ['varSummary_', varName];
    saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
    saveas(gcf,fullfile(folder,[saveFigName, '.png']));
    
end

%% Plot ptables
varNames = {'pMuTable', 'pMuTable_v', 'pMuTable_mu', 'pMuTable_sig', 'pSigTable', 'pSigTable_v'};
for k = 1:6
    
    
    varName = varNames{k};
    
    eval(['currentT = ', varName, ';']);
    
    maxSize = max(size(currentT));
    
    newT = zeros(maxSize);
    newT(1:end,1:end-1) = currentT;
    newT = newT + newT.';
    newT = newT + eye(size(newT))*0.5;
    
    eventNames = resultsTable.label;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    heatmap(eventNames, eventNames, newT);
    
    test = colormap(parula);
    test = ones(size(test));
    cutoff = 0.05;
    maxc = 0.5;
    minc = 0;
    cutoffIdx = ceil(size(test,1)*(cutoff)/(maxc - minc));
    test(1:cutoffIdx,2) = 0;
    test(1:cutoffIdx,3) = 0;
    test(1:cutoffIdx,1) = linspace(0.4,1,cutoffIdx);
    h = gca;
    
    colormap(test);
    title(varName);
    
    saveFigName = ['pVals_', varName];
    saveas(gcf,fullfile(folder,[saveFigName, '.fig']));
    saveas(gcf,fullfile(folder,[saveFigName, '.png']));
end



        