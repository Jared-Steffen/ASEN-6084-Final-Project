function [tracks,nodes] = TOMHT_AdaptiveQ(Imb,X,Y,P0,R,Q0,pD,rho_nu,rho_ex,T,alpha)
% TOMHT with compatible cluster hypotheses and score-threshold pruning
%
% Inputs:
%   Imb     : image sequence, only used for number of frames
%   X, Y    : cell arrays of detections per frame
%   P0      : initial covariance for births
%   R       : cell array, R{k}(:,:,i) = measurement covariance for meas i at frame k
%   Q       : 2x2 acceleration covariance for CV model
%   pD      : probability of detection
%   rho_nu  : new-target density
%   rho_ex  : extraneous measurement density
%   T       : pruning threshold
%
% Outputs:
%   tracks  : track structs
%   nodes   : node structs (hypothesis tree)

% ============================================================
% Initialize
% ============================================================
n = size(Imb,3);

nodes = struct( ...
    'id', {}, ...
    'track_id', {}, ...
    'parent_node_id', {}, ...
    'children_node_ids', {}, ...
    'frame', {}, ...
    'Z', {}, ...
    'P', {}, ...
    'Q', {}, ...
    'score', {}, ...
    'is_leaf', {}, ...
    'num_hits', {}, ...
    'num_misses', {}, ...
    'age', {}, ...
    'status', {});

tracks = struct( ...
    'id', {}, ...
    'root_node_id', {}, ...
    'leaf_node_ids', {});

clusters = struct( ...
    'id', {}, ...
    'track_ids', {}, ...
    'leaf_node_ids', {}, ...
    'measurement_ids', {}, ...            % cell array aligned with track_ids
    'unassigned_measurement_ids', {});    % vector of measurements that can spawn births

next_node_id = 1;
next_track_id = 1;
next_cluster_id = 1;

% CV pixel-space model
deltat = 1;
Phi = [eye(2), deltat.*eye(2);
       zeros(2), eye(2)];

Qkp1k = [deltat^3/3.*Q0, deltat^2/2.*Q0;
         deltat^2/2.*Q0, deltat.*Q0];

H = [eye(2), zeros(2)];

% ============================================================
% Main loop
% ============================================================
for k = 1:n

    Xk = cell2mat(X(k));
    Yk = cell2mat(Y(k));
    Rk = cell2mat(R(k));

    % --------------------------------------------------------
    % Reset cluster measurement bookkeeping for this frame
    % --------------------------------------------------------
    for c = 1:length(clusters)
        nT = length(clusters(c).track_ids);
        clusters(c).measurement_ids = cell(nT,1);
        clusters(c).unassigned_measurement_ids = [];
    end

    % --------------------------------------------------------
    % Frame 1 births
    % --------------------------------------------------------
    if k == 1
        for i = 1:length(Xk)

            Zk = [Xk(i); Yk(i); 0; 0];

            nodes(next_node_id).id = next_node_id;
            nodes(next_node_id).track_id = next_track_id;
            nodes(next_node_id).parent_node_id = 0;
            nodes(next_node_id).children_node_ids = [];
            nodes(next_node_id).frame = k;
            nodes(next_node_id).Z = Zk;
            nodes(next_node_id).P = P0;
            nodes(next_node_id).Q = Qkp1k;
            nodes(next_node_id).score = log(rho_nu / rho_ex);
            nodes(next_node_id).is_leaf = true;
            nodes(next_node_id).num_hits = 1;
            nodes(next_node_id).num_misses = 0;
            nodes(next_node_id).age = 1;
            nodes(next_node_id).status = "initialized";

            tracks(next_track_id).id = next_track_id;
            tracks(next_track_id).root_node_id = next_node_id;
            tracks(next_track_id).leaf_node_ids = next_node_id;

            clusters(next_cluster_id).id = next_cluster_id;
            clusters(next_cluster_id).track_ids = next_track_id;
            clusters(next_cluster_id).leaf_node_ids = {next_node_id};
            clusters(next_cluster_id).measurement_ids = {i};
            clusters(next_cluster_id).unassigned_measurement_ids = [];

            next_node_id = next_node_id + 1;
            next_track_id = next_track_id + 1;
            next_cluster_id = next_cluster_id + 1;
        end

        fprintf('Finished processing frame %d / %d\n', k, n);
        continue
    end

    % --------------------------------------------------------
    % Step 1: predict each leaf
    % --------------------------------------------------------
    predictions = struct( ...
        'node_id', {}, ...
        'track_id', {}, ...
        'Zhatm', {}, ...
        'Pm', {});

    p = 1;
    leaf_id_list = [];

    for i = 1:length(tracks)
        leaf_ids = tracks(i).leaf_node_ids;
        leaf_id_list = [leaf_id_list, leaf_ids];

        for j = 1:length(leaf_ids)
            leaf_id = leaf_ids(j);
            leaf_idx = find([nodes.id] == leaf_id, 1, 'first');
            if isempty(leaf_idx)
                continue
            end

            Zbarkm1 = nodes(leaf_idx).Z;
            Pbarkm1 = nodes(leaf_idx).P;
            Qhat = nodes(leaf_idx).Q;

            Zhatm = Phi * Zbarkm1;
            Pm = Phi * Pbarkm1 * Phi' + Qhat;

            predictions(p).node_id = leaf_id;
            predictions(p).track_id = tracks(i).id;
            predictions(p).Zhatm = Zhatm;
            predictions(p).Pm = Pm;
            p = p + 1;
        end
    end

    % --------------------------------------------------------
    % Step 2: gating
    % --------------------------------------------------------
    A = zeros(length(Xk), 1 + length(leaf_id_list));
    A(:,1) = 1;   % unassociated / birth column

    Pgate = 0.99;
    gamma = chi2inv(Pgate, 2);

    for i = 1:length(Xk)
        for j = 1:length(predictions)
            Zhatm = predictions(j).Zhatm;
            Pm = predictions(j).Pm;
            track_id = predictions(j).track_id;

            v = [Xk(i); Yk(i)] - H*Zhatm;
            S = H*Pm*H' + Rk(:,:,i);
            dM2 = v' * (S \ v);

            if dM2 < gamma
                A(i,j+1) = 1;

                % append measurement i to whatever cluster currently owns this track
                for m = 1:length(clusters)
                    loc = find(clusters(m).track_ids == track_id, 1, 'first');
                    if ~isempty(loc)
                        clusters(m).measurement_ids{loc} = ...
                            unique([clusters(m).measurement_ids{loc}, i]);
                    end
                end
            end
        end
    end

    % --------------------------------------------------------
    % Step 3: create measurement-only clusters for detections
    % that gated to nobody
    % --------------------------------------------------------
    new_clusters = struct( ...
        'id', {}, ...
        'track_ids', {}, ...
        'leaf_node_ids', {}, ...
        'measurement_ids', {}, ...
        'unassigned_measurement_ids', {});

    new_clusters_list = find(sum(A,2) == 1);

    for i = 1:length(new_clusters_list)
        meas_id = new_clusters_list(i);

        new_clusters(i).id = next_cluster_id;
        new_clusters(i).track_ids = [];
        new_clusters(i).leaf_node_ids = {};
        new_clusters(i).measurement_ids = {};
        new_clusters(i).unassigned_measurement_ids = meas_id;

        next_cluster_id = next_cluster_id + 1;
    end

    if ~isempty(new_clusters)
        clusters = [clusters, new_clusters];
    end

    % --------------------------------------------------------
    % Step 4: merge clusters that share tracks or measurements
    % --------------------------------------------------------
    nC = length(clusters);
    Adj = false(nC,nC);

    for i = 1:nC
        Adj(i,i) = true;
        for j = i+1:nC
            meas_i = cluster_measurement_union(clusters(i));
            meas_j = cluster_measurement_union(clusters(j));

            share_meas = ~isempty(intersect(meas_i, meas_j));
            share_track = ~isempty(intersect(clusters(i).track_ids, clusters(j).track_ids));

            if share_meas || share_track
                Adj(i,j) = true;
                Adj(j,i) = true;
            end
        end
    end

    Gc = graph(Adj);
    comp = conncomp(Gc);

    merged_clusters = struct( ...
        'id', {}, ...
        'track_ids', {}, ...
        'leaf_node_ids', {}, ...
        'measurement_ids', {}, ...
        'unassigned_measurement_ids', {});

    comps = unique(comp);
    mc = 1;

    for cc = 1:length(comps)
        idx = find(comp == comps(cc));

        merged_clusters(mc).id = mc;

        merged_track_ids = unique([clusters(idx).track_ids]);
        merged_clusters(mc).track_ids = merged_track_ids;

        nT = length(merged_track_ids);
        merged_leafs = cell(nT,1);
        merged_meas = cell(nT,1);

        for t = 1:nT
            tid = merged_track_ids(t);
            leaf_list = [];
            meas_list = [];

            for kk = idx
                loc = find(clusters(kk).track_ids == tid, 1, 'first');
                if ~isempty(loc)
                    leaf_list = [leaf_list, clusters(kk).leaf_node_ids{loc}];
                    meas_list = [meas_list, clusters(kk).measurement_ids{loc}];
                end
            end

            merged_leafs{t} = unique(leaf_list);
            merged_meas{t} = unique(meas_list);
        end

        merged_clusters(mc).leaf_node_ids = merged_leafs;
        merged_clusters(mc).measurement_ids = merged_meas;

        ua = [];
        for kk = idx
            ua = [ua, clusters(kk).unassigned_measurement_ids];
        end
        merged_clusters(mc).unassigned_measurement_ids = unique(ua);

        mc = mc + 1;
    end

    clusters = merged_clusters;

    % --------------------------------------------------------
    % Step 5: enumerate track-branch rows
    % --------------------------------------------------------
    nC = length(clusters);
    track_hyp_en = enumerate_cluster_hypotheses(clusters);

    % Build row_data from track_hyp_en
    row_data = cell(nC,1);

    for i = 1:nC
        current_cluster_hyp = track_hyp_en{i};
        nRows = size(current_cluster_hyp,1);

        row_data{i}(nRows) = struct( ...
            'llr', -Inf, ...
            'score', -Inf, ...
            'Zhatp', [], ...
            'Qhat', [], ...
            'Pp', [], ...
            'is_valid', false);

        for j = 1:nRows
            track_id = current_cluster_hyp(j,1);
            parent_node_id = current_cluster_hyp(j,2);
            meas_id = current_cluster_hyp(j,3);

            row_data{i}(j).llr = -Inf;
            row_data{i}(j).score = -Inf;
            row_data{i}(j).Zhatp = [];
            row_data{i}(j).Pp = [];
            row_data{i}(j).Qhat = [];
            row_data{i}(j).is_valid = false;

            % Existing track + measurement update
            if track_id ~= 0 && meas_id ~= 0
                pred_list_idx = find([predictions.node_id] == parent_node_id, 1, 'first');
                parent_node_idx = find([nodes.id] == parent_node_id, 1, 'first');

                if isempty(pred_list_idx) || isempty(parent_node_idx)
                    continue
                end

                parent_score = nodes(parent_node_idx).score;
                parent_Q = nodes(parent_node_idx).Q;

                G = [Xk(meas_id); Yk(meas_id)];
                Rassociated = Rk(:,:,meas_id);

                Zhatm = predictions(pred_list_idx).Zhatm;
                Pm = predictions(pred_list_idx).Pm;

                S = H*Pm*H' + Rassociated;
                Ki = Pm*H'/S;
                rhatm = G - H*Zhatm;
                Zhatp = Zhatm + Ki*rhatm;
                Pp = (eye(4) - Ki*H)*Pm*(eye(4) - Ki*H)' + Ki*Rassociated*Ki';

                llr = -0.5*(rhatm'*(S\rhatm)) ...
                      -0.5*log(det(2*pi*S)) ...
                      + log(pD) ...
                      - log(rho_ex);

                row_data{i}(j).llr = llr;
                row_data{i}(j).score = parent_score + llr;
                row_data{i}(j).Zhatp = Zhatp;
                row_data{i}(j).Pp = Pp;
                row_data{i}(j).Qhat = alpha.*parent_Q + (1-alpha) * Ki * (rhatm * rhatm') * Ki';
                row_data{i}(j).is_valid = true;

            % Existing track + missed detection
            elseif track_id ~= 0 && meas_id == 0
                pred_list_idx = find([predictions.node_id] == parent_node_id, 1, 'first');
                parent_node_idx = find([nodes.id] == parent_node_id, 1, 'first');

                if isempty(pred_list_idx) || isempty(parent_node_idx)
                    continue
                end

                parent_score = nodes(parent_node_idx).score;
                parent_Q = nodes(parent_node_idx).Q;

                Zhatp = predictions(pred_list_idx).Zhatm;
                Pp = predictions(pred_list_idx).Pm;

                llr = log(1-pD);

                row_data{i}(j).llr = llr;
                row_data{i}(j).score = parent_score + llr;
                row_data{i}(j).Zhatp = Zhatp;
                row_data{i}(j).Pp = Pp;
                row_data{i}(j).Qhat = parent_Q;
                row_data{i}(j).is_valid = true;

            % Birth row
            elseif track_id == 0 && meas_id ~= 0
                Zhatp = [Xk(meas_id); Yk(meas_id); 0; 0];
                Pp = P0;

                llr = log(rho_nu / rho_ex);

                row_data{i}(j).llr = llr;
                row_data{i}(j).score = llr;
                row_data{i}(j).Zhatp = Zhatp;
                row_data{i}(j).Pp = Pp;
                row_data{i}(j).Qhat = Qkp1k;
                row_data{i}(j).is_valid = true;
            end
        end
    end

    cluster_hyps = build_compatible_cluster_hypotheses(track_hyp_en, row_data);

    % --------------------------------------------------------
    % Step 6: instantiate child branches from surviving
    % compatible cluster hypotheses
    % --------------------------------------------------------
    children = struct( ...
        'cluster_id', {}, ...
        'parent_node_id', {}, ...
        'measurement_id', {}, ...
        'track_id', {}, ...
        'Zhatp', {}, ...
        'Pp', {}, ...
        'Qhat', {}, ...
        'llr', {}, ...
        'score', {});
    
    u = 1;
    
    for i = 1:nC
        hyps_i = cluster_hyps{i};
    
        if isempty(hyps_i)
            continue
        end
    
        % Continuing / missed-detection branches
        for h = 1:length(hyps_i)
            assigns = hyps_i(h).assignments;
    
            for a = 1:length(assigns)
                children(u).cluster_id = clusters(i).id;
                children(u).parent_node_id = assigns(a).parent_node_id;
                children(u).measurement_id = assigns(a).measurement_id;
                children(u).track_id = assigns(a).track_id;
                children(u).Zhatp = assigns(a).Zhatp;
                children(u).Pp = assigns(a).Pp;
                children(u).Qhat = assigns(a).Qhat;
                children(u).llr = assigns(a).llr;
                children(u).score = assigns(a).score;
                u = u + 1;
            end
        end
    
        % Birth branches:
        meas_union = cluster_measurement_union(clusters(i));
    
        % remove measurements already used by at least one kept hypothesis
        used_any = [];
        for h = 1:length(hyps_i)
            used_any = [used_any, hyps_i(h).used_measurements];
        end
        used_any = unique(used_any);
    
        birth_meas = setdiff(meas_union, used_any);
    
        for m = 1:length(birth_meas)
            meas_id = birth_meas(m);
    
            children(u).cluster_id = clusters(i).id;
            children(u).parent_node_id = 0;
            children(u).measurement_id = meas_id;
            children(u).track_id = next_track_id;
            children(u).Zhatp = [Xk(meas_id); Yk(meas_id); 0; 0];
            children(u).Pp = P0;
            children(u).Qhat = Qkp1k;
            children(u).llr = log(rho_nu / rho_ex);
            children(u).score = children(u).llr;
    
            next_track_id = next_track_id + 1;
            u = u + 1;
        end
    end
    
    % Remove repeats
    if ~isempty(children)
        child_keys = zeros(length(children),4);
    
        for ii = 1:length(children)
            if children(ii).parent_node_id == 0
                child_keys(ii,:) = [children(ii).cluster_id, -1, 0, children(ii).measurement_id];
            else
                child_keys(ii,:) = [children(ii).cluster_id, ...
                                    children(ii).track_id, ...
                                    children(ii).parent_node_id, ...
                                    children(ii).measurement_id];
            end
        end
    
        [~, idx] = unique(child_keys, 'rows', 'stable');
        children = children(idx);
    end
    % --------------------------------------------------------
    % Step 7: turn children into nodes
    % --------------------------------------------------------
    new_tracks = struct( ...
        'id', {}, ...
        'root_node_id', {}, ...
        'leaf_node_ids', {});

    new_id = 1;

    % wipe leaf lists before rebuilding from children
    for i = 1:length(tracks)
        tracks(i).leaf_node_ids = [];
    end
    for i = 1:length(clusters)
        clusters(i).leaf_node_ids = cell(length(clusters(i).track_ids),1);
        clusters(i).measurement_ids = cell(length(clusters(i).track_ids),1);
    end

    for i = 1:length(children)
        child_track = children(i).track_id;

        corresponding_node_idx = [];
        if children(i).parent_node_id ~= 0
            corresponding_node_idx = find([nodes.id] == children(i).parent_node_id, 1, 'first');
            if isempty(corresponding_node_idx)
                continue
            end
        end

        new_node_idx = length(nodes) + 1;

        nodes(new_node_idx).id = next_node_id;
        nodes(new_node_idx).track_id = child_track;
        nodes(new_node_idx).parent_node_id = children(i).parent_node_id;
        nodes(new_node_idx).children_node_ids = [];
        nodes(new_node_idx).frame = k;
        nodes(new_node_idx).Z = children(i).Zhatp;
        nodes(new_node_idx).P = children(i).Pp;
        nodes(new_node_idx).Q = children(i).Qhat;
        nodes(new_node_idx).score = children(i).score;
        nodes(new_node_idx).is_leaf = true;

        if children(i).measurement_id == 0 && children(i).parent_node_id ~= 0
            nodes(new_node_idx).num_hits = nodes(corresponding_node_idx).num_hits;
            nodes(new_node_idx).num_misses = nodes(corresponding_node_idx).num_misses + 1;
            nodes(new_node_idx).age = nodes(corresponding_node_idx).age + 1;
            nodes(new_node_idx).status = "coasting";

        elseif children(i).measurement_id ~= 0 && children(i).parent_node_id ~= 0
            nodes(new_node_idx).num_hits = nodes(corresponding_node_idx).num_hits + 1;
            nodes(new_node_idx).num_misses = nodes(corresponding_node_idx).num_misses;
            nodes(new_node_idx).age = nodes(corresponding_node_idx).age + 1;

            if nodes(new_node_idx).num_hits >= 5
                nodes(new_node_idx).status = "confirmed";
            else
                nodes(new_node_idx).status = "tentative";
            end

        else
            nodes(new_node_idx).num_hits = 1;
            nodes(new_node_idx).num_misses = 0;
            nodes(new_node_idx).age = 1;
            nodes(new_node_idx).status = "initialized";
        end

        % update parent bookkeeping
        if children(i).parent_node_id ~= 0 && ~isempty(corresponding_node_idx)
            nodes(corresponding_node_idx).is_leaf = false;
            nodes(corresponding_node_idx).children_node_ids = ...
                [nodes(corresponding_node_idx).children_node_ids, next_node_id];
        end

        corresponding_cluster_idx = find([clusters.id] == children(i).cluster_id, 1, 'first');

        % existing track child
        if children(i).parent_node_id ~= 0
            corresponding_track_idx = find([tracks.id] == child_track, 1, 'first');

            if ~isempty(corresponding_track_idx)
                tracks(corresponding_track_idx).leaf_node_ids = ...
                    [tracks(corresponding_track_idx).leaf_node_ids, next_node_id];
            end

            if ~isempty(corresponding_cluster_idx)
                loc = find(clusters(corresponding_cluster_idx).track_ids == child_track, 1, 'first');
                if ~isempty(loc)
                    clusters(corresponding_cluster_idx).leaf_node_ids{loc} = ...
                        [clusters(corresponding_cluster_idx).leaf_node_ids{loc}, next_node_id];

                    if children(i).measurement_id ~= 0
                        clusters(corresponding_cluster_idx).measurement_ids{loc} = ...
                            unique([clusters(corresponding_cluster_idx).measurement_ids{loc}, children(i).measurement_id]);
                    end
                end
            end

        % birth child
        else
            new_tracks(new_id).id = children(i).track_id;
            new_tracks(new_id).root_node_id = next_node_id;
            new_tracks(new_id).leaf_node_ids = next_node_id;
            new_id = new_id + 1;

            if ~isempty(corresponding_cluster_idx)
                clusters(corresponding_cluster_idx).track_ids(end+1) = children(i).track_id;
                clusters(corresponding_cluster_idx).leaf_node_ids{end+1} = next_node_id;
                clusters(corresponding_cluster_idx).measurement_ids{end+1} = children(i).measurement_id;
                clusters(corresponding_cluster_idx).unassigned_measurement_ids = ...
                    setdiff(clusters(corresponding_cluster_idx).unassigned_measurement_ids, children(i).measurement_id);
            end
        end

        next_node_id = next_node_id + 1;
    end

    if ~isempty(new_tracks)
        tracks = [tracks, new_tracks];
    end

    % --------------------------------------------------------
    % Step 8: Score-threshold pruning
    % --------------------------------------------------------
    if ~isempty(nodes)

        prune_mask = false(1, length(nodes));

        for ii = 1:length(nodes)
            if nodes(ii).age >= 5 && nodes(ii).score < T
                prune_mask(ii) = true;
            end
        end

        prune_ids = [nodes(prune_mask).id];

        % Remove pruned nodes from nodes
        if ~isempty(prune_ids)
            keep_mask = ~ismember([nodes.id], prune_ids);
            nodes = nodes(keep_mask);
        end

        if ~isempty(nodes)
            surviving_ids = [nodes.id];
        else
            surviving_ids = [];
        end

        % Rebuild children_node_ids
        for ii = 1:length(nodes)
            if ~isempty(nodes(ii).children_node_ids)
                nodes(ii).children_node_ids = ...
                    nodes(ii).children_node_ids( ...
                        ismember(nodes(ii).children_node_ids, surviving_ids) );
            end
        end

        % Recompute leaf flags
        if ~isempty(nodes)
            parent_ids = [nodes.parent_node_id];
            for ii = 1:length(nodes)
                nodes(ii).is_leaf = ~any(parent_ids == nodes(ii).id);
            end
        end

        % Update tracks.leaf_node_ids
        keep_track_mask = true(1, length(tracks));

        for ii = 1:length(tracks)
            if ~isempty(tracks(ii).leaf_node_ids)
                tracks(ii).leaf_node_ids = ...
                    tracks(ii).leaf_node_ids( ...
                        ismember(tracks(ii).leaf_node_ids, surviving_ids) );
            end

            if tracks(ii).root_node_id ~= 0 && ~ismember(tracks(ii).root_node_id, surviving_ids)
                keep_track_mask(ii) = false;
            elseif isempty(tracks(ii).leaf_node_ids)
                keep_track_mask(ii) = false;
            end
        end

        tracks = tracks(keep_track_mask);

        % Remove nodes belonging to deleted tracks
        if ~isempty(tracks)
            surviving_track_ids = [tracks.id];
            node_keep_mask = ismember([nodes.track_id], surviving_track_ids);
            nodes = nodes(node_keep_mask);
        else
            nodes = nodes([]);
            surviving_track_ids = [];
        end

        % Final cleanup after track-based node removal
        if ~isempty(nodes)
            surviving_ids = [nodes.id];
        else
            surviving_ids = [];
        end

        for ii = 1:length(nodes)
            if ~isempty(nodes(ii).children_node_ids)
                nodes(ii).children_node_ids = ...
                    nodes(ii).children_node_ids( ...
                        ismember(nodes(ii).children_node_ids, surviving_ids) );
            end
        end

        if ~isempty(nodes)
            parent_ids = [nodes.parent_node_id];
            for ii = 1:length(nodes)
                nodes(ii).is_leaf = ~any(parent_ids == nodes(ii).id);
            end
        end

        % Refresh track leaf lists
        for ii = 1:length(tracks)
            tid = tracks(ii).id;
            track_leaf_mask = [nodes.track_id] == tid & [nodes.is_leaf];
            tracks(ii).leaf_node_ids = [nodes(track_leaf_mask).id];
        end

        % Refresh clusters to match surviving tracks
        keep_cluster_mask = true(1, length(clusters));

        for ii = 1:length(clusters)
            keep_ct_mask = ismember(clusters(ii).track_ids, surviving_track_ids);

            clusters(ii).track_ids = clusters(ii).track_ids(keep_ct_mask);

            if ~isempty(clusters(ii).leaf_node_ids)
                clusters(ii).leaf_node_ids = clusters(ii).leaf_node_ids(keep_ct_mask);
            else
                clusters(ii).leaf_node_ids = {};
            end

            if ~isempty(clusters(ii).measurement_ids)
                clusters(ii).measurement_ids = clusters(ii).measurement_ids(keep_ct_mask);
            else
                clusters(ii).measurement_ids = {};
            end

            for jj = 1:length(clusters(ii).track_ids)
                tid = clusters(ii).track_ids(jj);
                cluster_leaf_mask = [nodes.track_id] == tid & [nodes.is_leaf];
                clusters(ii).leaf_node_ids{jj} = [nodes(cluster_leaf_mask).id];
            end

            % after pruning, measurement bookkeeping will be rebuilt next frame
            clusters(ii).measurement_ids = cell(length(clusters(ii).track_ids),1);
            clusters(ii).unassigned_measurement_ids = [];

            if isempty(clusters(ii).track_ids)
                keep_cluster_mask(ii) = false;
            end
        end

        clusters = clusters(keep_cluster_mask);
    end

    fprintf('Finished processing frame %d / %d\n', k, n);
end

end

% ============================================================
% Helpers
% ============================================================
function meas_union = cluster_measurement_union(cluster_i)
    meas_union = cluster_i.unassigned_measurement_ids;
    for q = 1:length(cluster_i.measurement_ids)
        meas_union = [meas_union, cluster_i.measurement_ids{q}];
    end
    meas_union = unique(meas_union);
end

function track_hyp_en = enumerate_cluster_hypotheses(clusters)
% For each cluster, return rows of:
%   [track_id, parent_node_id, measurement_id]
%
% Includes:
%   - one row per (track, leaf, measurement)
%   - one missed-detection row per (track, leaf)
%   - one birth row per measurement: [0, 0, meas_id]

    nC = length(clusters);
    track_hyp_en = cell(nC,1);

    for i = 1:nC
        cluster_tracks = clusters(i).track_ids;
        cluster_leaves = clusters(i).leaf_node_ids;
        meas_union = cluster_measurement_union(clusters(i));

        hyp_rows = zeros(0,3);

        for j = 1:length(cluster_tracks)
            track_id = cluster_tracks(j);
            leaf_ids = cluster_leaves{j};
            track_meas = [0, unique(meas_union)];   % 0 = miss

            for l = 1:length(leaf_ids)
                parent_node_id = leaf_ids(l);

                for m = 1:length(track_meas)
                    hyp_rows(end+1,:) = [track_id, parent_node_id, track_meas(m)];
                end
            end
        end

        % Birth rows: one per measurement in cluster
        for m = 1:length(meas_union)
            hyp_rows(end+1,:) = [0, 0, meas_union(m)];
        end

        track_hyp_en{i} = hyp_rows;
    end
end

function cluster_hyps = build_compatible_cluster_hypotheses(track_hyp_en, row_data)
% Inputs:
%   track_hyp_en : cell array, one cell per cluster
%                  each cell is Nr x 3 matrix:
%                  [track_id, parent_node_id, measurement_id]
%
%   row_data     : cell array, one cell per cluster
%                  row_data{i}(r) corresponds to track_hyp_en{i}(r,:)
%
% Output:
%   cluster_hyps : cell array, one cell per cluster
%                  each cell contains a struct array with fields:
%                       .assignments
%                       .score
%                       .used_measurements
%                       .unused_measurements

    nC = length(track_hyp_en);
    cluster_hyps = cell(nC,1);

    for i = 1:nC
        rows = track_hyp_en{i};

        if isempty(rows)
            cluster_hyps{i} = struct( ...
                'assignments', {}, ...
                'score', {}, ...
                'used_measurements', {}, ...
                'unused_measurements', {});
            continue
        end

        data_i = row_data{i};

        cont_mask = rows(:,1) ~= 0;
        cont_rows = rows(cont_mask,:);
        cont_data = data_i(cont_mask);

        if ~isempty(cont_data)
            valid_mask = [cont_data.is_valid];
            cont_rows = cont_rows(valid_mask,:);
            cont_data = cont_data(valid_mask);
        end

        cluster_meas = unique(rows(rows(:,3) ~= 0, 3));

        if isempty(cont_rows)
            cluster_hyps{i}(1).assignments = struct( ...
                'track_id', {}, ...
                'parent_node_id', {}, ...
                'measurement_id', {}, ...
                'llr', {}, ...
                'score', {}, ...
                'Zhatp', {}, ...
                'Pp', {}, ...
                'Qhat', {});
            cluster_hyps{i}(1).score = 0;
            cluster_hyps{i}(1).used_measurements = [];
            cluster_hyps{i}(1).unused_measurements = cluster_meas;
            continue
        end

        track_ids = unique(cont_rows(:,1));
        nT = length(track_ids);
        track_row_idxs = cell(nT,1);

        for t = 1:nT
            track_row_idxs{t} = find(cont_rows(:,1) == track_ids(t));
        end

        partial = struct( ...
            'row_idxs', {}, ...
            'used_meas', {}, ...
            'score', {});

        partial(1).row_idxs = [];
        partial(1).used_meas = [];
        partial(1).score = 0;

        for t = 1:nT
            candidate_idxs = track_row_idxs{t};

            new_partial = struct( ...
                'row_idxs', {}, ...
                'used_meas', {}, ...
                'score', {});

            u = 1;

            for h = 1:length(partial)
                for q = 1:length(candidate_idxs)
                    ridx = candidate_idxs(q);
                    meas_id = cont_rows(ridx,3);

                    if meas_id ~= 0 && any(partial(h).used_meas == meas_id)
                        continue
                    end

                    new_partial(u).row_idxs = [partial(h).row_idxs, ridx];

                    if meas_id == 0
                        new_partial(u).used_meas = partial(h).used_meas;
                    else
                        new_partial(u).used_meas = [partial(h).used_meas, meas_id];
                    end

                    new_partial(u).score = partial(h).score + cont_data(ridx).score;
                    u = u + 1;
                end
            end

            if isempty(new_partial)
                partial = struct( ...
                    'row_idxs', {}, ...
                    'used_meas', {}, ...
                    'score', {});
                break
            end

            % Keep only best X partial hypotheses
            [~,ord] = sort([new_partial.score], 'descend');
            keepN = min(8, length(ord));
            partial = new_partial(ord(1:keepN));
        
        end

        if isempty(partial)
            cluster_hyps{i} = struct( ...
                'assignments', {}, ...
                'score', {}, ...
                'used_measurements', {}, ...
                'unused_measurements', {});
            continue
        end

        hyps_i = struct( ...
            'assignments', {}, ...
            'score', {}, ...
            'used_measurements', {}, ...
            'unused_measurements', {});

        for h = 1:length(partial)
            ridxs = partial(h).row_idxs;

            assigns = struct( ...
                'track_id', {}, ...
                'parent_node_id', {}, ...
                'measurement_id', {}, ...
                'llr', {}, ...
                'score', {}, ...
                'Zhatp', {}, ...
                'Pp', {}, ...
                'Qhat', {});

            for a = 1:length(ridxs)
                rr = ridxs(a);

                assigns(a).track_id = cont_rows(rr,1);
                assigns(a).parent_node_id = cont_rows(rr,2);
                assigns(a).measurement_id = cont_rows(rr,3);
                assigns(a).llr = cont_data(rr).llr;
                assigns(a).score = cont_data(rr).score;
                assigns(a).Zhatp = cont_data(rr).Zhatp;
                assigns(a).Pp = cont_data(rr).Pp;
                assigns(a).Qhat = cont_data(rr).Qhat;
            end

            hyps_i(h).assignments = assigns;
            hyps_i(h).score = partial(h).score;
            hyps_i(h).used_measurements = unique(partial(h).used_meas);
            hyps_i(h).unused_measurements = setdiff(cluster_meas, partial(h).used_meas);
        end

        cluster_hyps{i} = hyps_i;
    end
end