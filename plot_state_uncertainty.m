function plot_state_uncertainty(t,P,sigma_scale)

% Formatting
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesGridAlpha',0.3)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultLineMarkerSize',8)


% 3-sigma bounds from covariance
N = length(t);
sig = zeros(N,4);
for k = 1:N
    sig(k,:) = sqrt(diag(P(:,:,k)))';
end
b = sigma_scale*sig;

figure();
labels = {'x', 'y', 'v_x', 'v_y'};
for i = 1:4
    subplot(4,1,i)
    line = plot(t,+b(:,i), '--', 'LineWidth', 3.0);
    hold on
    plot(t,-b(:,i), '--', 'LineWidth', 3.0,'Color',line.Color);
    grid on;
    xlabel('Frame');
    if i < 3
        ylim([-2.5 2.5])
    end
    ylabel(sprintf('Uncertainty in %s', labels{i}));
    title(sprintf('Pixel-Space Uncertainty in %s \\pm%g\\sigma', labels{i}, sigma_scale));
end


end