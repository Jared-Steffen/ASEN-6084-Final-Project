function [xhatp,Pp,Qhat,rhatm,rhatp,det_list] = LKF_AdaptiveQ(x0,P0,Imb,X,Y,R,Q0,alpha)
%{
Inputs:
    >xbar0: initailized state deviation
    >Pbar0: initialized state covariance matrix
    >Imb: Background subtracted image sequence
    >X: x detections for every frame
    >Y: y detections for every frame
    >R: measurement uncertainty matrix
    >Q: process noise covariance matrix
Outputs:
    >Xhat: full predicted state (dxhat + Xnom)
    >P: state history of estimated state covariance
    >rhatm: pre-fit measurement residuals
    >rhatp: post-fit residuals
%}


% Initialize
N = size(Imb,3);
xhatpkm1 = x0;
Ppkm1 = P0;
H = [eye(2) zeros(2)];
Phi = [eye(2) eye(2);
       zeros(2) eye(2)];
Gamma = [0.5.*eye(2);eye(2)];

% Iterative algorithm
for i = 1:N

    % Get STM
    if i > 1
        deltat = 1;
    else
        deltat = 0;
    end

    % Time update
    xhatm = Phi*xhatpkm1;

    % Do not add process noise before k =5 to gather needed information
    if i == 1
        Pm = Phi*Ppkm1*Phi + Gamma*Q0*Gamma';
    else
        Pm = Phi*Ppkm1*Phi' + Qhat(:,:,i-1);
    end

    % Get set of detections (x,y)
    Xi = cell2mat(X(i));
    Yi = cell2mat(Y(i));
    Ri = cell2mat(R(i));

    % Gate probability and chi-square threshold
    Pgate = 0.99;
    gamma = chi2inv(Pgate,2);

    % Find if any detections are within the gate
    dx = Xi - xhatm(1);
    dy = Yi - xhatm(2);
    dets = [];
    dM2_list = [];
    
    for k = 1:length(Xi)
        v = [Xi(k); Yi(k)] - H*xhatm;
        S = H*Pm*H' + Ri(:,:,k);
        dM2 = v' * (S \ v);
    
        if dM2 <= gamma
            dets(end+1) = k;
            dM2_list(end+1) = dM2;
        end
    end


    % If no associated detection, next time step
    if isempty(dets)
        det_list(i) = 0;
        xhatp(:,i) = xhatm;
        Pp(:,:,i) = Pm;
        rhatm(:,i) = NaN(2,1);
        rhatp(:,i) = NaN(2,1);
        xhatpkm1 = xhatm;
        Ppkm1 = Pm;
        Qhat(:,:,i) = Qhat(:,:,i-1);
        continue
    % If only one valid detection, associate
    elseif length(dets) == 1
        det_list(i) = 1;
        Xassociated = Xi(dets);
        Yassociated = Yi(dets);
        Rassociated = Ri(:,:,dets);
    % Otherwise, associate with smallest dM (Mahalanobis)
    else
        det_list(i) = length(dets);
        [~,idx] = min(dM2_list);
        det_min = dets(idx);
        Xassociated = Xi(det_min);
        Yassociated = Yi(det_min);
        Rassociated = Ri(:,:,det_min);
    end
    G = [Xassociated;Yassociated];

    % Measurement correction
    Si = H*Pm*H' + Rassociated;
    Ki = Pm*H'/Si;
    rhatm(:,i) = G - H*xhatm;
    xhatp(:,i) = xhatm + Ki*rhatm(:,i);
    Pp(:,:,i) = (eye(4) - Ki*H)*Pm*(eye(4) - Ki*H)' + Ki*Rassociated*Ki';
    rhatp(:,i) = G - H*xhatp(:,i);

    % Estimate process noise
    if i == 1
        Qhat(:,:,i) = Ki * (rhatm(:,i) * rhatm(:,i)') * Ki';
    else
        Qhat(:,:,i) = alpha * Qhat(:,:,i-1) + (1-alpha) * Ki * (rhatm(:,i) * rhatm(:,i)') * Ki';
    end

    % Move iteration forward
    Ppkm1 = Pp(:,:,i);
    xhatpkm1 = xhatp(:,i);
    
end

end
