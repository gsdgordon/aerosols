clc;
clear variables
close all;

folder = 'C:\Users\george\OneDrive - The University of Nottingham\SAVE\20201028\20202810_04_N60';
prefix = '--6-Run';


obscamTime_AT = datetime(2020,10,28,9,38,16);
AT_time = datetime(2020,10,28,8,47,54);
AT_diff = AT_time - obscamTime_AT;

obscamTime_N60 = datetime(2020,10,28,9,37,35);
N60_time = datetime(2020,10,28,8,48,56);
N60_diff = N60_time - obscamTime_N60;

N60_diff_AT = N60_diff - AT_diff;

[N60times, N60counts, N60diameters, N60shapes] = loadN60data(folder, prefix);