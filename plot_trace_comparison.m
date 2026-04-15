function plot_trace_comparison(P4, P3, p4_labels)

%{
Inputs:
    P4        : 3x3xN x M covariance matrix history
    P3        : 3x3xN covariance matrix history
    p4_labels : 1xM or Mx1 vector of values for legend entries
%}

[~,~,N,M] = size(P4);

trace4_pos = zeros(N,M);
trace3_pos = zeros(1,N);

trace4_vel = zeros(N,M);
trace3_vel = zeros(1,N);

% Compute traces for P4
for m = 1:M
    for i = 1:N
        trace4_pos(i,m) = trace(P4(1:2,1:2,i,m));
        trace4_vel(i,m) = trace(P4(3:4,3:4,i,m));
    end
end

% Compute trace for P3
for i = 1:N
    trace3_pos(i) = trace(P3(1:2,1:2,i));
    trace3_vel(i) = trace(P3(3:4,3:4,i));
end

%% Position
figure;
hold on;

h = gobjects(M+1,1);

for m = 1:M
    h(m) = plot(1:N, trace4_pos(:,m), 'o-', 'LineWidth', 1.5);
end

h(M+1) = plot(1:N, trace3_pos, 'o-', 'LineWidth', 1.5);

legend_entries = cell(M+1,1);
for m = 1:M
    legend_entries{m} = num2str(p4_labels(m));
end
legend_entries{M+1} = 'No Adaptive Q';

xlabel('Frame');
ylabel('Trace of Covariance');
title('Position Trace Comparison');
legend(h, legend_entries, 'Location', 'best');
grid on;

%% Velocity
figure;
hold on;

h = gobjects(M+1,1);

for m = 1:M
    h(m) = plot(1:N, trace4_vel(:,m), 'o-', 'LineWidth', 1.5);
end

h(M+1) = plot(1:N, trace3_vel, 'o-', 'LineWidth', 1.5);

legend_entries = cell(M+1,1);
for m = 1:M
    legend_entries{m} = num2str(p4_labels(m));
end
legend_entries{M+1} = 'No Adaptive Q';

xlabel('Frame');
ylabel('Trace of Covariance');
title('Velocity Trace Comparison');
legend(h, legend_entries, 'Location', 'best');
grid on;

end