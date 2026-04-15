clc; clear; close all

%% Matched Filter detection with Binary Hypothesis Testing

% Load in background subtracted image and other data
data = load("ImbData.mat");
Imb = data.Imb;
W = data.W;
sigma_n = data.sigma;
sigmaB = data.sigmaB;
sigma_psf = data.sigma_psf;

data2 = load("ImgSeq.mat");
ImgSeq = data2.ImgSeq;

% Get window size
h = floor(W/2);
m = W^2;

% Build PSF template
psf = zeros(W,W);
for r = 1:W
    for c = 1:W
        dy = r - (h+1);
        dx = c - (h+1);
        psf(r,c) = exp(-(dx^2 + dy^2)/(2*sigma_psf^2));
    end
end

% Get shat/gbar*
sbar = psf(:);
gstar = sbar / norm(sbar);

% Probability of false alarm
pFA = logspace(-6,-4,3);

% Centroid error constant values
var3 = 1/(12*m);

% Total signal variance 
varS = m.*sigma_n.^2;

%% MF algorithm
N = size(Imb,3);
for j = 1:length(pFA)
    for i = 1:N
        % Get threshold score
        TpFA(i,j) = norminv(1-pFA(j)).*(sqrt(sigma_n(i)^2+sigmaB(i)^2));

        % Get current image
        Imbi = Imb(:,:,i); 

        % Compute every pixel MF score
        g = reshape(gstar, W, W);
        MFscore  = conv2(Imbi, -g, 'same');

        % Threshold mask
        mask = (MFscore > TpFA(i,j));
        num_mask_pixels(i,j) = nnz(mask);

        % Initialize clustering outputs for this frame/pFA
        total_clusters(i,j) = 0;
        r = [];
        c = [];
        cluster_ids = [];
        valid_cluster_ids = [];

        % Decide if anything passed the test
        det = any(mask(:));

        % Cluster if an object is detected
        if det
            [r, c] = find(mask);
            X = [c r]; % indexes of detections
            cluster_ids = dbscan(X,1.5,3);
            valid_cluster_ids = unique(cluster_ids(cluster_ids ~= -1)); % exclude noise label
            total_clusters(i,j) = numel(valid_cluster_ids);
        end

        % Compute centroid estimates/error and pd for each detection
        xhatc_k = [];
        yhatc_k = [];
        Shat_k = [];
        varc_k = [];
        
        counter = 1;
        for kk = 1:numel(valid_cluster_ids)

            cid = valid_cluster_ids(kk);

            % Find pixel locations for each cluster
            rk = r(cluster_ids == cid);
            ck = c(cluster_ids == cid);

            % Place window on maximum MF score of cluster
            [max_MFscore,max_idx] = max(MFscore(sub2ind(size(MFscore),rk,ck)));
            r1 = max(1, rk(max_idx)-h);
            r2 = min(size(Imbi,1), rk(max_idx)+h);
            c1 = max(1, ck(max_idx)-h);
            c2 = min(size(Imbi,2), ck(max_idx)+h);

            window = Imbi(r1:r2, c1:c2);
            cols = c1:c2;
            rows = r1:r2;
            in_region1 = any(rows >= 1650 & rows <= 1920) && ...
                         any(cols >= 0 & cols < 190);

            in_region2 = any(rows >= 1500 & rows <= 1920) && ...
                         any(cols >= 900 & cols < 1080);

            if in_region1 || in_region2
                continue
            end

            den = sum(window,'all');
            xhatc_k(counter) = sum(cols .* sum(window,1)) / den;
            yhatc_k(counter) = sum(rows' .* sum(window,2)) / den;

            % Get additive noise error variance and total centroid variance
            Shat_k(kk) = sum(window,"all");
            var2 = sigma_n(i)^2/(Shat_k(counter)^2 + m*sigma_n(i)^2)*(W*sum((-h:h).^2,"all"));
            varc_k(counter) = var2 + var3;
            
            counter = counter + 1;
        end

        % Save for current frame,pFa combo
        xhatc{i,j} = xhatc_k;
        yhatc{i,j} = yhatc_k;
        Shat{i,j} = Shat_k;
        varc{i,j} = varc_k;

    end
end

%% Plots

figure();
subplot(1,2,1)
plot(1:N,num_mask_pixels,'o-')
xlabel('Frame Number')
ylabel('Number of Pixel Detections')
title('Number of Individual Detections for Various Probability of False Alarms')
legend(compose("pFA = %.3g", pFA))
subplot(1,2,2)
plot(1:N,total_clusters,'o-')
xlabel('Frame Number')
ylabel('Number of Clusters')
title('Number of Clusters for Various Probability of False Alarms')
legend(compose("pFA = %.3g", pFA))

% for i = 1:size(Imb,3)
%     figure();
    % subplot(1,3,1)
    % imagesc(Imb(:,:,i))
    % hold on
    % colormap gray
    % plot(xhatc{i,1},yhatc{i,1},'.','MarkerSize',13)
    % subplot(1,3,2)
    % imagesc(Imb(:,:,i))
    % hold on
    % colormap gray
    % plot(xhatc{i,6},yhatc{i,6},'.','MarkerSize',13)
    % title('Estimated Cluster Centroids for Various Probability of False Alarms for Frame 1')
    % subplot(1,3,3)
    % imshow(Imb(:,:,i))
    % hold on
    % colormap gray
    % plot(xhatc{i,3},yhatc{i,3},'.','MarkerSize',13)
% end


figure();
subplot(1,3,1)
imagesc(ImgSeq(:,:,1))
hold on
colormap gray
plot(xhatc{1,1},yhatc{1,1},'.','MarkerSize',13)
subplot(1,3,2)
imagesc(ImgSeq(:,:,1))
hold on
colormap gray
plot(xhatc{1,2},yhatc{1,2},'.','MarkerSize',13)
title('Estimated Cluster Centroids for Various Probability of False Alarms for Frame 1')
subplot(1,3,3)
imagesc(ImgSeq(:,:,1))
hold on
colormap gray
plot(xhatc{1,3},yhatc{1,3},'.','MarkerSize',13)

%----------------------------------------------
figure();
subplot(1,3,1)
imshow(Imb(:,:,50))
hold on
colormap gray
plot(xhatc{50,1},yhatc{50,1},'.','MarkerSize',13)
subplot(1,3,2)
imshow(Imb(:,:,50))
hold on
colormap gray
plot(xhatc{50,2},yhatc{50,2},'.','MarkerSize',13)
title('Estimated Cluster Centroids for Various Probability of False Alarms for Frame 1')
subplot(1,3,3)
imshow(Imb(:,:,50))
hold on
colormap gray
plot(xhatc{50,3},yhatc{50,3},'.','MarkerSize',13)

%----------------------------------------------
figure();
subplot(1,3,1)
imshow(Imb(:,:,end))
hold on
colormap gray
plot(xhatc{end,1},yhatc{end,1},'.','MarkerSize',13)
subplot(1,3,2)
imshow(Imb(:,:,end))
hold on
colormap gray
plot(xhatc{end,2},yhatc{end,2},'.','MarkerSize',13)
title('Estimated Cluster Centroids for Various Probability of False Alarms for Frame 1')
subplot(1,3,3)
imshow(Imb(:,:,end))
hold on
colormap gray
plot(xhatc{end,3},yhatc{end,3},'.','MarkerSize',13)

save('DetectionData2.mat','xhatc','yhatc','varc','pFA')