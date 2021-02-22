% Aggegrate N60 data
clc;
clear variables;
close all;

diameters = [];
shapes = [];

N60_area = 12606*7181; %0.58 magnification, in micronrs
N60_depth = 50000;
N60_vol = N60_area/(1e6^2) * N60_depth/(1e6);

for k=1:4
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k), '.mat']);
    diameters = cat(1,diameters,N60diameters);
    shapes = cat(1,shapes,N60shapes);
end

diamBins = logspace(log10(10),log10(400),10);
avDiam = sqrt((diamBins(1:end-1)).*(diamBins(2:end)));
shapeBins = 0:0.1:1;

diameters = diameters(diameters ~= 0);
shapes = shapes(shapes ~= 0);

param = lognfit(diameters);

subplot(2,1,1)
counts = histcounts(diameters, diamBins);
counts = counts/N60_vol/4;
histogram('BinCounts', counts, 'BinEdges', diamBins);
xlabel('diameters \mum');
ylabel('no. particles /m^3 per procedure')
set(gca, 'XScale', 'log')


subplot(2,1,2)
counts = histcounts(shapes, shapeBins);
counts = counts/N60_vol/4;
histogram('BinCounts', counts, 'BinEdges', shapeBins);
xlabel('shape factors');
ylabel('no. particles /m^3 per procedure')


for k=[1,2,4]
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_corr_', num2str(k), '.mat']);
    corrCurve = circshift(corrCurve, -3);
    if k == 1
        sumCorr = corrCurve;
    else
        sumCorr = sumCorr + corrCurve;
    end
end

sumCorr = sumCorr/3;

figure;
plot(lags, sumCorr);
xlabel('delay in seconds');
ylabel('normalised correlation');
title('Correlation of N60 with AeroTrak 10um');

diameters_pos = [];
shapes_pos = [];
for k=1:5
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k) '_Position changes.mat']);
    corrCurve = circshift(corrCurve, -3);
    diameters_pos = cat(1,diameters_pos,N60diameters);
    shapes_pos = cat(1,shapes_pos,N60shapes);
end
diameters_pos = diameters_pos(diameters_pos ~= 0);
shapes_pos = shapes_pos(shapes_pos ~= 0);

diameters_D2 = [];
shapes_D2 = [];
for k=1:4
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k) '_D2 reached.mat']);
    corrCurve = circshift(corrCurve, -3);
    diameters_D2 = cat(1,diameters_D2,N60diameters);
    shapes_D2 = cat(1,shapes_D2,N60shapes);
end
diameters_D2 = diameters_D2(diameters_D2 ~= 0);
shapes_D2 = shapes_D2(shapes_D2 ~= 0);

diameters_fun = [];
shapes_fun = [];
for k=1:4
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k) '_Fundal retroflexion.mat']);
    corrCurve = circshift(corrCurve, -3);
    diameters_fun = cat(1,diameters_fun,N60diameters);
    shapes_fun = cat(1,shapes_fun,N60shapes);
end
diameters_fun = diameters_fun(diameters_fun ~= 0);
shapes_fun = shapes_fun(shapes_fun ~= 0);


diameters_ext = [];
shapes_ext = [];
for k=1:4
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k) '_Procedure ends.mat']);
    corrCurve = circshift(corrCurve, -3);
    diameters_ext = cat(1,diameters_ext,N60diameters);
    shapes_ext = cat(1,shapes_ext,N60shapes);
end
diameters_ext = diameters_ext(diameters_ext ~= 0);
shapes_ext = shapes_ext(shapes_ext ~= 0);

diameters_cou = [];
shapes_cou = [];
for k=[2,4]
    load(['C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\N60_agg_', num2str(k) '_Cough.mat']);
    corrCurve = circshift(corrCurve, -3);
    diameters_cou = cat(1,diameters_cou,N60diameters);
    shapes_cou = cat(1,shapes_cou,N60shapes);
end
diameters_cou = diameters_cou(diameters_cou ~= 0);
shapes_cou = shapes_cou(shapes_cou ~= 0);


figure;
subplot(2,3,1)
counts = histcounts(diameters_pos, diamBins);
counts = counts/N60_vol/5/60;
histogram('BinCounts', counts, 'BinEdges', diamBins);
xlabel('diameters \mum');
ylabel('no particles per m^3 / s');
set(gca, 'XScale', 'log')
title('Position changes');
ylim([0,6e3]);
ylim1 = ylim;


subplot(2,3,2)
counts = histcounts(diameters_D2, diamBins);
counts = counts/N60_vol/4/30;
histogram('BinCounts', counts, 'BinEdges', diamBins);
xlabel('diameters \mum');
ylabel('no particles per m^3 / s');
set(gca, 'XScale', 'log')
title('D2 reached');
ylim(ylim1);

subplot(2,3,3)
counts = histcounts(diameters_fun, diamBins);
counts = counts/N60_vol/4/60;
histogram('BinCounts', counts, 'BinEdges', diamBins);
xlabel('diameters \mum');
ylabel('no particles per m^3 / s');
set(gca, 'XScale', 'log')
title('Fundal retroflexion');
ylim(ylim1);

subplot(2,3,4)
counts = histcounts(diameters_ext, diamBins);
counts = counts/N60_vol/4/60;
histogram('BinCounts', counts, 'BinEdges', diamBins);
xlabel('diameters \mum');
ylabel('no particles per m^3 / s');
set(gca, 'XScale', 'log')
title('Final extubation');
ylim(ylim1);

subplot(2,3,5)
histogram(diameters_cou, diamBins);
xlabel('diameters \mum');
set(gca, 'XScale', 'log')
title('Coughing');
ylim(ylim1);

%% Shapes
figure;
subplot(2,3,1)
counts = histcounts(shapes_pos, shapeBins);
counts = counts/N60_vol/5/60;
histogram('BinCounts', counts, 'BinEdges', shapeBins);
xlabel('shape factor');
ylabel('no particles per m^3 / s');
title('Position changes');
ylim([0,6e3]);
ylim1 = ylim;


subplot(2,3,2)
counts = histcounts(shapes_D2, shapeBins);
counts = counts/N60_vol/4/60;
histogram('BinCounts', counts, 'BinEdges', shapeBins);
xlabel('shape factor');
ylabel('no particles per m^3 / s');
title('D2 reached');
ylim(ylim1);

subplot(2,3,3)
counts = histcounts(shapes_fun, shapeBins);
counts = counts/N60_vol/5/60;
histogram('BinCounts', counts, 'BinEdges', shapeBins);
xlabel('shape factor');
ylabel('no particles per m^3 / s');
title('Fundal retroflexion');
ylim(ylim1);

subplot(2,3,4)
counts = histcounts(shapes_ext, shapeBins);
counts = counts/N60_vol/5/60;
histogram('BinCounts', counts, 'BinEdges', shapeBins);
xlabel('shape factor');
ylabel('no particles per m^3 / s');
title('Final extubation');
ylim(ylim1);

subplot(2,3,5)
histogram(shapes_cou, shapeBins);
xlabel('shape factor');
title('Coughing');
ylim(ylim1);
ylabel('# per m^3');


