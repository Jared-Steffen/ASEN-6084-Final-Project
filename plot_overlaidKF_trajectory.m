function plot_overlaidKF_trajectory(t,Imb,Xhat,detections)

% Get initial and final images
Imb1 = Imb(:,:,1);

figure();
subplot(1,2,1)
imagesc(Imb1)
hold on
colormap gray
plot(Xhat(1,:),Xhat(2,:),'.-','MarkerSize',13)
title('Overlaid Track (on 1st Image)')
subplot(1,2,2)
plot(t,detections)
xlabel('Frame')
ylabel('Number of Detections')
title('Number of Detections in Gate')

end

