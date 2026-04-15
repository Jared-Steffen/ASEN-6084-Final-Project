function final_tracks = extractBestTracks(tracks, nodes)
    
final_tracks = struct( ...
        'track_id', {}, ...
        'best_leaf_node_id', {}, ...
        'best_leaf_score', {}, ...
        'node_ids', {}, ...
        'frames', {}, ...
        'states', {}, ...
        'covariances', {}, ...
        'scores', {}, ...
        'num_hits', {}, ...
        'num_misses', {}, ...
        'ages', {}, ...
        'statuses', {});

    if isempty(tracks) || isempty(nodes)
        return
    end

    all_node_ids = [nodes.id];
    out_idx = 1;

    for i = 1:length(tracks)

        leaf_ids = tracks(i).leaf_node_ids;

        if isempty(leaf_ids)
            continue
        end

        % Only keep valid leaf IDs
        valid_leaf_mask = ismember(leaf_ids, all_node_ids);
        leaf_ids = leaf_ids(valid_leaf_mask);

        if isempty(leaf_ids)
            continue
        end

        % ------------------------------------------------------------
        % Find highest-score leaf
        % ------------------------------------------------------------
        leaf_node_indices = zeros(1,length(leaf_ids));
        leaf_scores = -inf(1,length(leaf_ids));

        for j = 1:length(leaf_ids)
            idx = find(all_node_ids == leaf_ids(j), 1, 'first');
            leaf_node_indices(j) = idx;
            leaf_scores(j) = nodes(idx).score;
        end

        [best_leaf_score, best_j] = max(leaf_scores);

        % 🔴 NEW: skip low-quality tracks
        if best_leaf_score <= 0
            continue
        end

        best_leaf_idx = leaf_node_indices(best_j);
        best_leaf_id = nodes(best_leaf_idx).id;

        % ------------------------------------------------------------
        % Trace path
        % ------------------------------------------------------------
        path_node_ids = [];
        path_frames = [];
        path_states = [];
        path_covariances = [];
        path_scores = [];
        path_num_hits = [];
        path_num_misses = [];
        path_ages = [];
        path_statuses = strings(0,1);

        current_idx = best_leaf_idx;

        while true
            current_node = nodes(current_idx);

            path_node_ids(end+1) = current_node.id;
            path_frames(end+1) = current_node.frame;
            path_states(:,end+1) = current_node.Z;
            path_covariances(:,:,size(path_covariances,3)+1) = current_node.P;
            path_scores(end+1) = current_node.score;
            path_num_hits(end+1) = current_node.num_hits;
            path_num_misses(end+1) = current_node.num_misses;
            path_ages(end+1) = current_node.age;
            path_statuses(end+1,1) = string(current_node.status);

            if current_node.parent_node_id == 0
                break
            end

            parent_idx = find(all_node_ids == current_node.parent_node_id, 1, 'first');
            if isempty(parent_idx)
                break
            end

            current_idx = parent_idx;
        end

        % Reverse order (root → leaf)
        path_node_ids = fliplr(path_node_ids);
        path_frames = fliplr(path_frames);
        path_states = fliplr(path_states);
        path_covariances = path_covariances(:,:,end:-1:1);
        path_scores = fliplr(path_scores);
        path_num_hits = fliplr(path_num_hits);
        path_num_misses = fliplr(path_num_misses);
        path_ages = fliplr(path_ages);
        path_statuses = flipud(path_statuses);

        % Store
        final_tracks(out_idx).track_id = tracks(i).id;
        final_tracks(out_idx).best_leaf_node_id = best_leaf_id;
        final_tracks(out_idx).best_leaf_score = best_leaf_score;
        final_tracks(out_idx).node_ids = path_node_ids;
        final_tracks(out_idx).frames = path_frames;
        final_tracks(out_idx).states = path_states;
        final_tracks(out_idx).covariances = path_covariances;
        final_tracks(out_idx).scores = path_scores;
        final_tracks(out_idx).num_hits = path_num_hits;
        final_tracks(out_idx).num_misses = path_num_misses;
        final_tracks(out_idx).ages = path_ages;
        final_tracks(out_idx).statuses = path_statuses;

        out_idx = out_idx + 1;
    end
end

