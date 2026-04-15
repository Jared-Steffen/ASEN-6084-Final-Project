clc; clear; close all

%% Load in imaging data

% Setup
scriptPath = fileparts(mfilename('fullpath'));
folder = fullfile(scriptPath, 'image_seq');
files = dir(fullfile(folder, '*.png'));
N = length(files);

% Load images
i = 1;
for k = 1:N
    filename = fullfile(folder, files(k).name);
    I = imread(filename); 
    I = rgb2gray(I);
    I = im2double(I);    
    ImgSeq(:,:,i) = I;
    i = i + 1;
end

save('ImgSeq.mat','ImgSeq','-v7.3')

% Input video
% vidFile = "IMG_1455_lossless.mp4";
% 
% v = VideoReader(vidFile);
% 
% % Read first frame to get dimensions
% v.CurrentTime = 0;
% frame = read(v,1);
% imshow(frame)
% gray = rgb2gray(frame);
% gray = im2double(gray);
% 
% [ny,nx] = size(gray);
% 
% % Estimate number of frames
% N = floor(v.Duration * v.FrameRate);
% ImgSeq = zeros(ny, nx, N);
% 
% % Store first frame
% ImgSeq(:,:,1) = gray;
% 
% k = 1;
% while hasFrame(v)
%     k = k + 1;
% 
%     frame = readFrame(v);
%     gray = im2double(rgb2gray(frame));
% 
%     ImgSeq(:,:,k) = gray;
% end
% 
% % Trim in case frame estimate was too large
% ImgSeq = ImgSeq(:,:,1:k);
% 
% % Save
% save('ImgSeq.mat','ImgSeq','-v7.3')
