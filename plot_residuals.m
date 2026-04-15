function plot_residuals(t,rhatm,rhatp)

% Remove NaNs for histogram/statistics use
rmx = rhatm(1,~isnan(rhatm(1,:)));
rmy = rhatm(2,~isnan(rhatm(2,:)));
rpx = rhatp(1,~isnan(rhatp(1,:)));
rpy = rhatp(2,~isnan(rhatp(2,:)));

figure;

% Time-history plots
subplot(4,2,1)
plot(t, rhatm(1,:), '.', 'MarkerSize', 13)
xlabel('Frame')
ylabel('Residual Value')
title('x Position Pre-Fit Residual')
grid on

subplot(4,2,2)
plot(t, rhatp(1,:), '.', 'MarkerSize', 13)
xlabel('Frame')
ylabel('Residual Value')
title('x Position Post-Fit Residual')
grid on

subplot(4,2,3)
plot(t, rhatm(2,:), '.', 'MarkerSize', 13)
xlabel('Frame')
ylabel('Residual Value')
title('y Position Pre-Fit Residual')
grid on

subplot(4,2,4)
plot(t, rhatp(2,:), '.', 'MarkerSize', 13)
xlabel('Frame')
ylabel('Residual Value')
title('y Position Post-Fit Residual')
grid on

% Histograms + Gaussian overlays
subplot(4,2,5)
plot_hist_with_gaussian(rmx)
title('x Position Pre-Fit Histogram')

subplot(4,2,6)
plot_hist_with_gaussian(rpx)
title('x Position Post-Fit Histogram')

subplot(4,2,7)
plot_hist_with_gaussian(rmy)
title('y Position Pre-Fit Histogram')

subplot(4,2,8)
plot_hist_with_gaussian(rpy)
title('y Position Post-Fit Histogram')

end


function plot_hist_with_gaussian(data)

histogram(data,20, 'Normalization', 'pdf');
hold on

mu = mean(data);
sigma = std(data);

if sigma > 0
    x = linspace(min(data), max(data), 200);
    y = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((x-mu)/sigma).^2);
    plot(x, y, 'LineWidth', 2)
end

xlabel('Residual Value')
ylabel('PDF')
grid on
hold off

end
