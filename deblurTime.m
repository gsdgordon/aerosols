clc;
clear variables;
close all;

T = imread('vlcsnap-2020-10-14-11h25m57s762.png');

T = double(T);
T = mean(T,3);
T_crop = T;%T(600:end,1600:end);

imagesc(T_crop);
axis image;
colormap gray;

p1 = [831, 627];
p2 = [802, 609];

vec = p1- p2;
len = norm(vec);
len = 5;
theta = -1*atan2(vec(2),vec(1))*180/pi;
theta = 95;
theta_test = theta;
len_test = len;
ISR = 0.0005;

figure;
figure('units','normalized','outerposition',[0 0 1 1])

%for len_test = 2:20
    PSF = fspecial('motion',(len_test),theta_test);
    imagesc(PSF);
    axis image;

    T_de = deconvwnr(T_crop,PSF, ISR);
    imagesc(T_de);
    colormap gray;
    axis image;
    axis off;
    text(20,20,['var = ' num2str(len_test)])
    
    pause(1);
%end

figure;
subplot(1,2,1)
imagesc(T_crop);
colormap gray;
axis image;
axis off;
title('raw image');

subplot(1,2,2)
imagesc(T_de);
colormap gray;
axis image;
axis off;
title('corrected image');


