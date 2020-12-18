clc;
clear variables;
close all;

testVals = linspace(-2,14,2000);

for k = 1:size(testVals,2)-1
    dx(k) = exp(testVals(k+1)) - exp(testVals(k));
end


plot(testVals(1:end-1), log(dx))