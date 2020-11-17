clc;
clear variables;
close all;

load('resultsTable.mat');

noiseMean = exp(resultsTable.mean_n(2));
noiseMean = 0;
noiseStd = exp(resultsTable.std_n(2));
idx1 = 5;
idx2 = 6;
temp_n = resultsTable.n_raw{idx1};
dataSet1 = temp_n(:);
temp_n = resultsTable.n_raw{idx2};
dataSet2 = temp_n(:);
computeSignificance(dataSet1, dataSet2, noiseMean, noiseStd, resultsTable.mean_n(idx1), resultsTable.std_n(idx1), resultsTable.mean_n(idx2), resultsTable.std_n(idx2));