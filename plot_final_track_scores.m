function plot_final_track_scores(final_tracks)
% plot_final_track_scores
%
% Plots score vs frame for each final track and visually denotes status.
%
% Expected fields in each final_tracks(i):
%   .track_id
%   .frames
%   .scores
%   .statuses
%
% Status styles:
%   initialized -> black square
%   tentative   -> blue circle
%   coasting    -> orange triangle
%   confirmed   -> green diamond

    if isempty(final_tracks)
        error('final_tracks is empty.');
    end

    nTracks = length(final_tracks);

    % Choose subplot layout
    nCols = ceil(sqrt(nTracks));
    nRows = ceil(nTracks / nCols);

    figure;
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Status display settings
    status_names   = ["initialized", "tentative", "coasting", "confirmed"];
    marker_styles  = {'s', 'o', '^', 'd'};
    marker_colors  = [ ...
        1.0 1.0 0.0;   % initialized - yellow
        0.0 0.4470 0.7410; % tentative - blue
        0.8500 0.3250 0.0980; % coasting - orange
        0.4660 0.6740 0.1880]; % confirmed - green

    for i = 1:nTracks
        nexttile;
        hold on;
        grid on;

        frames = final_tracks(i).frames(:)';
        scores = final_tracks(i).scores(:)';

        % Convert statuses to string array if needed
        statuses = final_tracks(i).statuses;
        if ischar(statuses)
            statuses = string(cellstr(statuses));
        elseif iscell(statuses)
            statuses = string(statuses);
        else
            statuses = string(statuses);
        end

        % Base line through all scores
        plot(frames, scores, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.25);

        % Overlay markers by status
        for s = 1:length(status_names)
            mask = strcmpi(strtrim(statuses), status_names(s));
            if any(mask)
                plot(frames(mask), scores(mask), ...
                    marker_styles{s}, ...
                    'MarkerSize', 6, ...
                    'LineWidth', 1.2, ...
                    'MarkerEdgeColor', marker_colors(s,:), ...
                    'MarkerFaceColor', marker_colors(s,:));
            end
        end

        xlabel('Frame');
        xlim([0 85])
        ylabel('Score');
        title(sprintf('Track %d', final_tracks(i).track_id));

        % Make legend once per subplot only for statuses that appear
        legend_entries = {};
        legend_handles = [];

        % Add dummy handles in desired order
        for s = 1:length(status_names)
            mask = strcmpi(strtrim(statuses), status_names(s));
            if any(mask)
                h = plot(nan, nan, marker_styles{s}, ...
                    'MarkerSize', 6, ...
                    'LineWidth', 1.2, ...
                    'MarkerEdgeColor', marker_colors(s,:), ...
                    'MarkerFaceColor', marker_colors(s,:));
                legend_handles(end+1) = h;
                legend_entries{end+1} = char(status_names(s)); 
            end
        end

        if ~isempty(legend_handles)
            legend(legend_handles, legend_entries, 'Location', 'best');
        end
    end
end
