function plot_cov_trace_comparison(final_tracks_1, final_tracks_2, track_map)
% plot_cov_trace_comparison
%
% Inputs:
%   final_tracks_1 : struct array of tracks from set 1
%   final_tracks_2 : struct array of tracks from set 2
%   track_map      : Nx2 matrix
%                    col 1 = track_id in final_tracks_1
%                    col 2 = matching track_id in final_tracks_2
%
% Plots, for each mapped pair, the trace of the covariance matrices
% versus frame on the same plot using o- markers.
    nPairs = size(track_map,1);

    %% Position
    figure;
    tiledlayout('flow');

    for i = 1:nPairs
        tid1 = track_map(i,1);
        tid2 = track_map(i,2);

        idx1 = find([final_tracks_1.track_id] == tid1, 1);
        idx2 = find([final_tracks_2.track_id] == tid2, 1);

        nexttile; hold on; grid on;

        % Normal
        if ~isempty(idx1)
            covs1 = final_tracks_1(idx1).covariances;
            frames1 = final_tracks_1(idx1).frames;

            n1 = size(covs1,3) - 1; % ignore last
            if n1 > 0
                vals = zeros(1,n1);
                for k = 1:n1
                    vals(k) = trace(covs1(1:2,1:2,k));
                end
                plot(frames1(1:n1), vals, 'o-', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Track %d', tid1));
            end
        end

        % Adaptive Q
        if ~isempty(idx2)
            covs2 = final_tracks_2(idx2).covariances;
            frames2 = final_tracks_2(idx2).frames;

            n2 = size(covs2,3) - 1;
            if n2 > 0
                vals = zeros(1,n2);
                for k = 1:n2
                    vals(k) = trace(covs2(1:2,1:2,k));
                end
                plot(frames2(1:n2), vals, '--o', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Adaptive Q Track %d', tid2));
            end
        end

        title(sprintf('Tracks %d vs %d', tid1, tid2));
        xlabel('Frame'); ylabel('P diag');
        legend('Location','best');
    end

    %% Velocity
    figure;
    tiledlayout('flow');

    for i = 1:nPairs
        tid1 = track_map(i,1);
        tid2 = track_map(i,2);

        idx1 = find([final_tracks_1.track_id] == tid1, 1);
        idx2 = find([final_tracks_2.track_id] == tid2, 1);

        nexttile; hold on; grid on;

        % Normal
        if ~isempty(idx1)
            covs1 = final_tracks_1(idx1).covariances;
            frames1 = final_tracks_1(idx1).frames;

            n1 = size(covs1,3) - 1;
            if n1 > 0 && size(covs1,1) >= 4
                vals = zeros(1,n1);
                for k = 1:n1
                    vals(k) = trace(covs1(3:4,3:4,k));
                end
                plot(frames1(1:n1), vals, 'o-', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Track %d', tid1));
            end
        end

        % Adaptive Q
        if ~isempty(idx2)
            covs2 = final_tracks_2(idx2).covariances;
            frames2 = final_tracks_2(idx2).frames;

            n2 = size(covs2,3) - 1;
            if n2 > 0 && size(covs2,1) >= 4
                vals = zeros(1,n2);
                for k = 1:n2
                    vals(k) = trace(covs2(3:4,3:4,k));
                end
                plot(frames2(1:n2), vals, '--o', 'LineWidth', 1.2, ...
                    'DisplayName', sprintf('Adaptive Q Track %d', tid2));
            end
        end

        title(sprintf('Tracks %d vs %d', tid1, tid2));
        xlabel('Frame'); ylabel('P diag');
        legend('Location','best');
    end
end