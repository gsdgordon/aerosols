clc;
clear variables;
close all;

load('resultsTable.mat');

noiseMean = exp(resultsTable.mean_v(2));
noiseMean = 0;
noiseStd = exp(resultsTable.std_v(2));
idx1 = 4;
idx2 = 5;
temp_n = resultsTable.v_raw{idx1};
dataSet1 = temp_n(:);
temp_n = resultsTable.v_raw{idx2};
dataSet2 = temp_n(:)+10e-12;
dataSet2 = lognrnd(resultsTable.mean_v(idx2), resultsTable.std_v(idx2), [100,1]) + normrnd(noiseMean, noiseStd,[100,1]);
computeSignificance(dataSet1, dataSet2, noiseMean, noiseStd, resultsTable.mean_v(idx1), resultsTable.std_v(idx1), resultsTable.mean_v(idx2), resultsTable.std_v(idx2), 'muMinIn', -34);