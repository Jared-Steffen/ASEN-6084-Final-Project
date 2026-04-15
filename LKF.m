function [xhatp,Pp,rhatm,rhatp,det_list] = LKF(x0,P0,Imb,X,Y,R,Q)
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

% Iterative algorithm
for i = 1:N

    % Get STM
    if i > 1
        deltat = 1;
    else
        deltat = 0;
    end
    Phi = [eye(2) deltat.*eye(2);
           zeros(2) eye(2)];

    % Get DT Process Noise 
    Qkp1k = [deltat^3/3.*Q deltat^2/2.*Q;
             deltat^2/2.*Q deltat.*Q];

    % Time update
    xhatm = Phi*xhatpkm1;
    Pm = Phi*Ppkm1*Phi' + Qkp1k;

    % Get set of detections (x,y)
    Xi = cell2mat(X(i));
    Yi = cell2mat(Y(i));
    Ri = cell2mat(R(i));

    % % Set gate 
    % r = 5*sqrt(Pm(1,1)+Pm(2,2));

    % Gate probability and chi-square threshold
    Pgate = 0.99;
    gamma = chi2inv(Pgate,2);

    % Find if any detections are within the gate
    dx = Xi - xhatm(1);
    dy = Yi - xhatm(2);
    % inside_gate = (dx.^2+dy.^2) < r^2;
    % dets = find(inside_gate);
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
        continue
    % If only one valid detection, associate
    elseif length(dets) == 1
        det_list(i) = 1;
        Xassociated = Xi(dets);
        Yassociated = Yi(dets);
        Rassociated = Ri(:,:,dets);
    % Otherwise, associate with smallest dM (Mahalanobis)
    else
        % det_list(i) = length(dets);
        % for k = 1:length(dets)
        %     det = dets(k);
        %     S = H*Pm*H' + Ri(:,:,det);
        %     v = [Xi(det); Yi(det)] - (H*xhatm);
        %     dM(k) = sqrt(v' * (S \ v));
        % end
        % [~,idx] = min(dM);
        % det_min = dets(idx);
        % Xassociated = Xi(det_min);
        % Yassociated = Yi(det_min);
        % Rassociated = Ri(:,:,det_min);
        det_list(i) = length(dets);
        [~,idx] = min(dM2_list);
        det_min = dets(idx);
        Xassociated = Xi(det_min);
        Yassociated = Yi(det_min);
        Rassociated = Ri(:,:,det_min);
    end
    G = [Xassociated;Yassociated];

    % Measurement correction
    Ki = Pm*H'/(H*Pm*H' + Rassociated);
    rhatm(:,i) = G - H*xhatm;
    xhatp(:,i) = xhatm + Ki*rhatm(:,i);
    Pp(:,:,i) = (eye(4) - Ki*H)*Pm*(eye(4) - Ki*H)' + Ki*Rassociated*Ki';
    rhatp(:,i) = G - H*xhatp(:,i);

    % Move iteration forward
    Ppkm1 = Pp(:,:,i);
    xhatpkm1 = xhatp(:,i);
    
end

end
