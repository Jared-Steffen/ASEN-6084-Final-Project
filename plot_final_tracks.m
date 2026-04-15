function plot_final_tracks(I, final_tracks, frame_idx, show_labels)

    if nargin < 3 || isempty(frame_idx)
        frame_idx = inf;
    end

    if nargin < 4 || isempty(show_labels)
        show_labels = true;
    end

    figure;
    imagesc(I);
    hold on;

    if ndims(I) == 2
        colormap gray;
    end

    title('Final Tracks Overlaid on Image');
    xlabel('x [pixels]');
    ylabel('y [pixels]');

    nTracks = length(final_tracks);
    if nTracks == 0
        hold off;
        return
    end

    cmap = lines(max(nTracks,7));

    for i = 1:nTracks

        frames_i = final_tracks(i).frames;
        states_i = final_tracks(i).states;

        if isempty(frames_i) || isempty(states_i)
            continue
        end

        keep = frames_i <= frame_idx;

        if ~any(keep)
            continue
        end

        states_plot = states_i(:,keep);

        x = states_plot(1,:);
        y = states_plot(2,:);

        c = cmap(i,:);

        hp = plot(x, y, '-', 'LineWidth', 1.5, 'Color', c);
        plot(x(1), y(1), 'o', 'MarkerSize', 6, 'LineWidth', 1.2, ...
            'Color', c, 'MarkerEdgeColor', c);
        plot(x(end), y(end), 'x', 'MarkerSize', 8, 'LineWidth', 1.5, ...
            'Color', c);

        if show_labels
            text(x(end)+2, y(end)+2, sprintf('%d', final_tracks(i).track_id), ...
                'Color', hp.Color, 'FontSize', 20, 'FontWeight', 'bold');
        end
    end

    hold off;
end