function [ pwProb ] = pt_compute_spatial_pairwise_probabilities( p, dets, cidxs, pwIdxsAllrel1, graph, ...
                                                      hist_pairwise, nextreg, spatial_model, ...
                                                      use_graph_subset)

un_weight = 1;
pw_weight = 1;

addEpsHist = 100;
negHistMin = 0;
posHistMin = 0;

num_proposals = size(dets.unProb, 1);
pwProb = nan(floor(0.5*num_proposals^2*length(cidxs))+num_proposals^2,5);

[q1,q2] = meshgrid(1:num_proposals,1:num_proposals);
idxsAllrel = [q1(:) q2(:)];
idxsAll = idxsAllrel;

idxsExc = idxsAll(:,1) >= idxsAll(:,2);
idxsSameAll = idxsAll(~idxsExc,:);
idxSameAllrel = idxsAllrel(~idxsExc,:);

pid1 = dets.partClass(idxsSameAll(:,1));
pid2 = dets.partClass(idxsSameAll(:,2));

idxs = (pid1 == pid2);
idxsSame = idxsSameAll(idxs,:);
idxsDiff = idxsSameAll(~idxs,:);

xs = dets.unPos(:, 1);
ys = dets.unPos(:, 2);
ss = dets.scale;

if hist_pairwise
    half_size = 33;
    boxes = [xs-half_size, ys-half_size, xs+half_size, ys+half_size];
end

same_part_pw = isfield(p, 'same_part_locref_pw') && p.same_part_locref_pw;

idxStart = 1;
if ~same_part_pw
    % take scale of either of these to compensate for scale variation.
    % Multiplying with scale should help
    % since detections at the same scale will have better binaries than
    % on different scales.
    s = ss(idxsSame(:,1)); 

    x1 = xs(idxsSame(:,1)) .* s;
    y1 = ys(idxsSame(:,1)) .* s;

    x2 = xs(idxsSame(:,2)) .* s;
    y2 = ys(idxsSame(:,2)) .* s;

    dists = sqrt((x1 - x2).^2 +(y1 - y2).^2);
    prob = 1./(1+exp(0.2*dists-7.5));
else
    error('This mode is not supported at the moment. FixMe: add checking with respect to cidxs');
%     locations = dets.unPos+squeeze(dets.locationRefine(:,k,:));
%     featSame = get_spatial_features_same_part_regr_img(locations,idxsSame);
%     featSame_augm = get_augm_spatial_features_same_regr(featSame);
%     feat_norm = getFeatNorm(featSame_augm,same.training_opts.X_min,same.training_opts.X_max); 
% 
%     ex = sparse(double(feat_norm));
%     model = spatial_model_same{cidx}.log_reg;
%     [~,acc,prob_out] = predict(zeros(size(ex,1),1), ex, model, '-b 1');
%     %fprintf('same part prob old neg %f pos %f\n', prob(4911), prob(2806));
%     prob = prob_out(:,1); 
    %if cidx == 14
    %    fprintf('same part prob new neg %f pos %f\n', prob(4911), prob(2806));
    %    assert(false);
    %end
end

idxs = idxStart:idxStart+size(idxsSame,1)-1;
pwProb(idxs,1) = dets.index(idxsSame(:,1));
pwProb(idxs,2) = dets.index(idxsSame(:,2));
pwProb(idxs,3:4) = repmat(dets.partClass(idxsSame(:,1)), 1,2);
pwProb(idxs,5) = prob(:,1);

idxStart = idxStart + size(idxsSame,1);

fprintf('compute pairwise across parts prob\n');
idx_pw_prob = 0;
num_pw_prob = 1;
pw_prob_out = zeros(size(idxsAllrel, 1), num_pw_prob);

for m = 1:length(cidxs)-1
    for n = m+1:length(cidxs)
        cidx1 = cidxs(m);
        cidx2 = cidxs(n);
        for sidx1 = 1:length(pwIdxsAllrel1)
            if (pwIdxsAllrel1{sidx1}(1)== cidx1 && pwIdxsAllrel1{sidx1}(2)== cidx2)
                break;
            end
        end
        assert(sidx1 <= length(pwIdxsAllrel1));
        
        idxs = bitand(dets.partClass(idxsDiff(:,1)) == cidx1, dets.partClass(idxsDiff(:,2)) == cidx2);
        idxsRel = idxsDiff(idxs,:);
        
        idxsOpp = bitand(dets.partClass(idxsDiff(:,1)) == cidx2, dets.partClass(idxsDiff(:,2)) == cidx1);
        idxsRelOpp = idxsDiff(idxsOpp,:);
        idxsRelOpp = idxsRelOpp(:,[2,1]);
        idxsRel = cat(1,idxsRel, idxsRelOpp);
        
        if(isempty(idxsRel))
            continue;
        end

        is_neighbor = false;%isfield(p, 'allpairs') && p.allpairs;

        if nextreg % && is_neighbor
            [is_neighbor, forward_edge, backward_edge] = search_in_graph(graph, cidx1, cidx2);
            if isfield(p, 'neighbor_locref') && p.neighbor_locref
                featDiff = get_spatial_features_neighbour_locref(dets.unPos,idxsRel, ...
                                                                 squeeze(nextReg(:,[forward_edge backward_edge],:)), ...
                                                                 squeeze(locationRefine(:,[cidx1 cidx2],:)));
                featDiff = get_augm_spatial_features_diff_neighbour_locref(featDiff);
            else
                featDiff = pt_get_spatial_features_neighbour_img(dets.unPos,idxsRel,squeeze(dets.nextReg(:,[forward_edge backward_edge],:)), dets.scale);
                featDiff_augm = get_augm_spatial_features_diff_neighbour(featDiff, p);
            end
            feat_norm = getFeatNorm(featDiff_augm,spatial_model{sidx1}.diff.training_opts.X_min,spatial_model{sidx1}.diff.training_opts.X_max); 

        elseif hist_pairwise
            rotOffset = p.rot_offset;
            [featDiff,~,cb1Diff,cb2Diff] = get_spatial_features_diff_img_hist_dist(boxes,idxsRel,rotOffset);

            negHist = spatial_model{sidx1}.diff.negHist;
            negHist = negHist + addEpsHist;
            negHist = negHist./(sum(sum(negHist)));
            negHist(negHist < negHistMin) = negHistMin;
            spatial_model{sidx1}.diff.negHist = negHist;

            posHist = spatial_model{sidx1}.diff.posHist;
            posHist = posHist + addEpsHist;
            posHist = posHist./(sum(sum(posHist)));
            posHist(posHist < posHistMin) = posHistMin;
            spatial_model{sidx1}.diff.posHist = posHist;

            prob_out = computePosterior2D(featDiff,spatial_model{sidx1}.diff.edges1,spatial_model{sidx1}.diff.edges2, ...
                spatial_model{sidx1}.diff.posHist,spatial_model{sidx1}.diff.negHist,...
                spatial_model{sidx1}.diff.nPos,spatial_model{sidx1}.diff.nNeg);

        else
            rotOffset = p.rot_offset;
            [featDiff,~,cb1Diff,cb2Diff] = get_spatial_features_diff_img_dx_dy_dense(unPos,idxsRel,rotOffset,dets.unProb);
            featDiff = get_augm_spatial_features_diff_dx_dy(featDiff, spatial_model{sidx1}.diff.training_opts.X_pos_mean, p);
            feat_norm = getFeatNorm(featDiff,spatial_model{sidx1}.diff.training_opts.X_min,spatial_model{sidx1}.diff.training_opts.X_max);
        end

        if ~hist_pairwise || (is_neighbor && nextreg)
            ex = sparse(double(feat_norm));
            model = spatial_model{sidx1}.diff.log_reg;
            if p.liblinear_predict
                [~,~,prob_out] = predict(zeros(size(ex,1),1), ex, model, '-b 1');
                prob = prob_out(:,1);
            else
                x = [ full(ex) ones(size(ex, 1), 1)];
                prob = 1./(1+exp(-x*model.w'));
            end
        end

        if use_graph_subset
            [~, idx_in_sub] = ismember([cidx1 cidx2], graph_subset, 'rows');
            if idx_in_sub == 0
                prob_out = 0.5;
            end
        end

        prob = prob.^pw_weight;

        idxs = idxStart:idxStart+size(idxsRel,1)-1;
        pwProb(idxs,1) = dets.index(idxsRel(:,1));
        pwProb(idxs,2) = dets.index(idxsRel(:,2));
        pwProb(idxs,3) = dets.partClass(idxsRel(:,1));
        pwProb(idxs,4) = dets.partClass(idxsRel(:,2));
        pwProb(idxs,5) = prob;

        if (cidx1 == 7 && cidx2 == 8 && false)
        %if is_neighbor && false
            %forward_edge, backward_edge
            % idx = find(ismember(idxsAllrel(:,1),[29,53,13]) & ismember(idxsAllrel(:,2),[19 32 33 34]));
            dispProb = pwProb(idxs,5);
            %prob_cidx1= unProb(idxsAll(:, 1), cidx1);
            %prob_cidx2= unProb(idxsAll(:, 2), cidx2);
            %dispProb = dispProb .* prob_cidx1 .* prob_cidx2;
            dispProb = dispProb;

            [val,idx] = sort(dispProb,'descend');
            num_to_show = 20;
            idxsOracle = idx(1:min(length(idx),num_to_show));
            val = pwProb(idxs(idxsOracle),5);
            val(1:num_to_show, :)
            figure(101);clf;
            img = imread(annolist(i).image.name);
            imagesc(img);
            hold on;
            for j = 1:length(idxsOracle)
                %axis equal;
                %cb1 = cb1Diff(idxsOracle(j),:);
                %cb2 = cb2Diff(idxsOracle(j),:);

                cb1 = unPos(idxsRel(idxsOracle(j),1), :);
                cb2 = unPos(idxsRel(idxsOracle(j),2), :);

                plot([cb1(1); cb2(1)],[cb1(2); cb2(2)],'g-','lineWidth',1);
                %line([cb1(1); cb2(1)],[cb1(2); cb2(2)],'Color','g','lineWidth',3);
                plot(cb1(1), cb1(2),'bo','MarkerSize',5,'MarkerFaceColor','b');
                plot(cb2(1), cb2(2),'ro','MarkerSize',5,'MarkerFaceColor','r');
%                             [~,ind1pos] = histc(featDiff(idxsOracle(j),1),spatial_model{sidx1}.diff.edges1);
%                             [~,ind2pos] = histc(featDiff(idxsOracle(j),2),spatial_model{sidx1}.diff.edges2);
%                             [~,ind1neg] = histc(featDiff(idxsOracle(j),1),spatial_model{sidx1}.diff.edges1);
%                             [~,ind2neg] = histc(featDiff(idxsOracle(j),2),spatial_model{sidx1}.diff.edges2);
                %fprintf('prob: %1.2f, dist: %1.2f, rot: %1.2f, cidx1: %1.2f, cidx2: %1.2f\n',val(j),...
                %    featDiff(idx(j),1),featDiff(idx(j),2)*180/pi,featDiff(idx(j),offsetIdx+cidx1),featDiff(idx(j),offsetIdx+length(cidxs)+cidx2));
            end
            figure(11);
            histogram(prob);
            pause;
        end

        idxStart = idxStart + size(idxsRel,1);
    end
end

%return;

pwProb(idxStart-1:end,:) = [];
%pwProb(pwProb(:,5) <= 0.2, 5) = 0;
end

