clc; clear; close all

% Load data
data = load("DetectionData2.mat");
pFA = data.pFA;
varc = data.varc;
xhatc = data.xhatc;
yhatc = data.yhatc;

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

% Process Noise Covariance
Q = diag([1 1]);

% Initial covariance for each track
P0 = diag([4 4 25 25]);

% New object and clutter statistics
lambda_c = pFA(pFA_level)*length(reshape(ImgSeq(:,:,1),[],1));
lambda_nu = 3/85;
lambda_ex = lambda_nu+lambda_c;

V = length(reshape(ImgSeq(:,:,1),[],1)); 
rho_ex = (lambda_c + lambda_nu) / V;
rho_nu = lambda_nu/V;

% Track probability of detection
pD = 0.95;

[tracks, nodes] = TOMHT(ImgSeq,X,Y,P0,R,Q,pD,rho_nu,rho_ex,0);

% Extract tracks
final_tracks = extractBestTracks(tracks, nodes);
plot_final_tracks(ImgSeq(:,:,1), final_tracks, [], true)
plot_final_tracks(ImgSeq(:,:,end), final_tracks, [], true)
plot_final_track_scores(final_tracks)

%% Adaptive Q
[tracks_adQ,nodes_adQ] = TOMHT_AdaptiveQ(ImgSeq,X,Y,P0,R,Q,pD,rho_nu,rho_ex,0,0.8);

% Extract tracks
final_tracks_adQ = extractBestTracks(tracks_adQ, nodes_adQ);
plot_final_tracks(ImgSeq(:,:,1), final_tracks_adQ, [], true)
plot_final_tracks(ImgSeq(:,:,end), final_tracks_adQ, [], true)
plot_final_track_scores(final_tracks_adQ)

%% Compare
track_map = [1,1;3,3;4,4;5,5;9,9;15,18;36,33];
plot_cov_trace_comparison(final_tracks, final_tracks_adQ, track_map)
