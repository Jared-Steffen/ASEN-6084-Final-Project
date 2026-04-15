clc; clear; close all

%% Single Target Tracking

% Load Data
data = load('DetectionData2.mat');
xhatc = data.xhatc;
yhatc = data.yhatc;
varc = data.varc;
pFA = data.pFA;

data2 = load("ImgSeq.mat");
ImgSeq = data2.ImgSeq;

% Select data
pFA_level = 3;
X = xhatc(:,pFA_level);
Y = yhatc(:,pFA_level);
varC = varc(:,pFA_level);

% Form measurment uncertainty matrices
for i = 1:length(varC)
    varCi = varC{i};
    for j = 1:length(varCi)
        Rclust(:,:,j) = diag([varCi(j), varCi(j)]);
    end
    R{i} = Rclust;
end


% Initialize KF
% x0 = [85,1495,0,0]';
x0 = [416,1890,0,0]';
P0 = diag([4 4 25 25]);
Q = diag([1 1]);
Gamma = [0.5.*eye(2);eye(2)];
Qk = Gamma*Q*Gamma';

% KF Run
[xhatp,Pp,rhatm,rhatp,det_list] = LKF(x0,P0,ImgSeq,X,Y,R,Q);

%% Plots
N =  size(ImgSeq,3);

plot_state_uncertainty(1:N,Pp,3)
plot_overlaidKF_trajectory(1:N,ImgSeq,xhatp,det_list)
plot_residuals(1:N,rhatm,rhatp)


alpha = 0.5:0.1:0.9;
for i = 1:length(alpha)
    [xhatp_adQ(:,:,i),Pp_adQ(:,:,:,i),Qhat(:,:,:,i),rhatm(:,:,i),rhatp(:,:,i),det_list(:,:,i)] = LKF_AdaptiveQ(x0,P0,ImgSeq,X,Y,R,Q,alpha(i));
    plot_state_uncertainty(1:N,Pp_adQ(:,:,:,i),3)
    plot_overlaidKF_trajectory(1:N,ImgSeq,xhatp_adQ(:,:,i),det_list(:,:,i))
    plot_residuals(1:N,rhatm(:,:,i),rhatp(:,:,i))
end

plot_trace_comparison(Pp_adQ, Pp,alpha)

figure();
subplot(1,2,1)
imagesc(ImgSeq(:,:,1))
hold on
colormap gray
plot(xhatp(1,:),xhatp(2,:),'.-','MarkerSize',13)
plot(xhatp_adQ(1,:,4),xhatp_adQ(2,:,4),'.-','MarkerSize',13)
title('Overlaid Track Comparison (on 1st Image)')
subplot(1,2,2)
imagesc(ImgSeq(:,:,end))
hold on
colormap gray
plot(xhatp(1,:),xhatp(2,:),'.-','MarkerSize',13)
plot(xhatp_adQ(1,:,4),xhatp_adQ(2,:,4),'.-','MarkerSize',13)
title('Overlaid Track Comparison (on Final Image)')