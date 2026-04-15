clc; clear; close all

%% Subtract Background for each frame

% PSF width (6*sigma_psf) -- information from iPhone metadata
lambda = 550*1e-9; % wavelngth of light (assumption) [nm->m]
p = 1.22*1e-6; % pixel pitch [micro-m->m]
N = 1.78; % focal number
f = 6.765*1e-3; % aperature length [mm->m]
theta = 1.22*lambda*N; % Airy-disk small angle approximation
sigma_psf = theta/(3*p);
W = ceil(6*sigma_psf)+1;

% Load image sequence
data = load("ImgSeq.mat");
ImgSeq = data.ImgSeq;

% Subtract background
[Imb,Bhat,sigma,sigmaB] = BackgroundSubtraction(ImgSeq,sigma_psf);

%% Plots
figure();
subplot(1,3,1)
imshow(ImgSeq(:,:,1))
grid off
title('Original Image')
subplot(1,3,2)
imshow(Bhat(:,:,1))
grid off
title('Estimated Background')
subplot(1,3,3)
imshow(Imb(:,:,1))
grid off
title('Background Subtracted Image')

figure();
subplot(1,3,1)
imshow(ImgSeq(:,:,50))
grid off
title('Original Image')
subplot(1,3,2)
imshow(Bhat(:,:,50))
grid off
title('Estimated Background')
subplot(1,3,3)
imshow(Imb(:,:,50))
grid off
title('Background Subtracted Image')


figure();
subplot(1,3,1)
imshow(ImgSeq(:,:,end))
grid off
title('Original Image')
subplot(1,3,2)
imshow(Bhat(:,:,end))
grid off
title('Estimated Background')
subplot(1,3,3)
imshow(Imb(:,:,end))
grid off
title('Background Subtracted Image')

% Other plot
M = size(ImgSeq,3);

figure();
subplot(1,2,1)
plot(1:M,sigmaB)
xlabel('Frame Number')
ylabel('Background Estimate Uncertainty')
title('Background Estimate Uncertainty for Image Sequence')
subplot(1,2,2)
plot(1:M,sigma)
xlabel('Frame Number')
ylabel('Read Noise Uncertainty')
title('Read Noise Uncertainty for Image Sequence')

save('ImbData.mat','Imb','Bhat','sigma','sigmaB','sigma_psf','W')


