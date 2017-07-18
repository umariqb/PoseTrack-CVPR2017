function [X_pos, X_neg] = pt_get_spatial_features_dm(expidx,cidx)

RandStream.setGlobalStream ...
        (RandStream('mt19937ar','seed',42));

p = pt_exp_params(expidx);
bVis = false;
dist_thresh = 5;
pidxs = p.pidxs;
[~,parts] = util_get_parts24();
nFeatSample = p.nFeatSample;

save_file = [p.ptPairwiseDir sprintf('/feat_spatial_idx_%d.mat', cidx)];

if exist(save_file, 'file') == 2
    fprintf('Loading feature file: %s\n', save_file);
    load(save_file);
    fprintf('Loaded %s\n',save_file);
    return;
end

fprintf('Loading annotations from: %s\n', p.trainGT);
load(p.trainGT, 'annolist');



num_videos = length(annolist);

num_frames = 0;
num_persons = 0;
for s=1:num_videos
    num_frames  = num_frames  + annolist(s).num_frames;
    num_persons = num_persons + annolist(s).num_persons;        
end

dim = 5;
% initialize the arrays, otherwise too slow when using cat
X_pos = zeros(num_frames*num_persons*nFeatSample,dim,'single');
X_neg = zeros(num_frames*num_persons*nFeatSample,dim,'single');

% example counter for each class
n_pos = 0;
n_neg = 0;
for vIdx=1:num_videos
    fprintf('Computing temporal features.\t Video %d/%d\n', vIdx, num_videos);
    ann = annolist(vIdx);
    vid_name    = ann.name;
    vid_dir = fullfile(p.vidDir, vid_name);
    num_frames  = ann.num_frames;
    num_persons = ann.num_persons;        
    fn = dir([vid_dir,'/*.jpg']);
    
    corres_dir = fullfile(p.correspondences, vid_name);
    scoremaps_dir = fullfile(p.scoremaps, vid_name);
    
    if num_persons == 1
        nFeatSamplePers = nFeatSample;
    else
        nFeatSamplePers = uint16(nFeatSample/5);
    end

    for fIdx=1:num_frames
            
        fr_fn = fullfile(vid_dir, fn(fIdx).name);
        [~,fr_name,~] = fileparts(fr_fn);

        corres_fn = fullfile(corres_dir, [fr_name1,'_',fr_name2,'.txt']);
        corres_pts = pt_load_dm_correspondences(corres_fn);
        corres_pts1 = corres_pts(:,1:2);
        corres_pts2 = corres_pts(:,3:4);
        
        score_fn1 = fullfile(scoremaps_dir, fr_name1);
        score_fn2 = fullfile(scoremaps_dir, fr_name2);
        
        if(p.multiscale || true)
            locations1 = pt_generate_location_candidates(score_fn1, p, cidx);
            locations2 = pt_generate_location_candidates(score_fn2, p, cidx);
        end
        
        for pIdx = 1:num_persons
            annopoints1 = ann.annopoints{pIdx, pair(1)};
            if(isempty(annopoints1))
                continue;
            end
            joints = pt_get_anno_joints(ann.annopoints{pIdx, pair(1)}, pidxs, parts);
            gt_joint1 = joints(cidx,:);
            
            for ppIdx = 1:num_persons
                
                annopoints2 = ann.annopoints{ppIdx,pair(2)};
                
                if(isempty(annopoints2))
                    continue;
                end
                
                joints = pt_get_anno_joints(ann.annopoints{ppIdx,pair(2)}, pidxs, parts);
                gt_joint2 = joints(cidx,:);
                
                if (isnan(gt_joint1(1)) || isnan(gt_joint2(1)))
                    continue;
                end
                
                dx = locations1(:,1)- gt_joint1(1);
                dy = locations1(:,2)- gt_joint1(2);
                dists_to_gt1 = sqrt(dx.^2 + dy.^2);

                dx = locations2(:,1)- gt_joint2(1);
                dy = locations2(:,2)- gt_joint2(2);
                dists_to_gt2 = sqrt(dx.^2 + dy.^2);
                
                % compute dist_thresh according to scales
                scales = locations1(:,3);
                dist_thresh_all1 = scales * dist_thresh;
                scales = locations2(:,3);
                dist_thresh_all2 = scales * dist_thresh;

                idxs_bbox1 = find(dists_to_gt1 <= dist_thresh_all1);
                idxs_bbox2 = find(dists_to_gt2 <= dist_thresh_all2);
                
                % perform nms to combine detections from different scales
                % the one with better score will be kept
                locs1 = locations1(idxs_bbox1, :);
                idxs_bbox1 = nms_distance(locations1(idxs_bbox1, :), 3);
                
                locs2 = locations2(idxs_bbox2, :);         
                idxs_bbox2 = nms_distance(locations2(idxs_bbox2, :), 3);
                
                [i,j] = meshgrid(idxs_bbox1, idxs_bbox2);
                idxAll = [i(:) j(:)];

                if(pIdx == ppIdx) % positive features
                    nFeat = min(nFeatSamplePers,size(idxAll,1));
                    idxs_rnd = randperm(size(idxAll,1));
                    idxsPair = idxs_rnd(1:nFeat);
                    idxs_bbox_pair = idxAll(idxsPair,:);  
                    feat = pt_get_temporal_features_img_dm(p, locs1, corres_pts1, locs2, corres_pts2, idxs_bbox_pair);
                    idxs = n_pos+1:n_pos+size(feat,1);
                    X_pos(idxs,:) = feat;
                    n_pos = n_pos + size(feat,1);
                else  % negative features
                    nFeat = uint16(min(nFeatSamplePers/3,size(idxAll,1)));
                    idxs_rnd = randperm(size(idxAll,1));
                    idxsPair = idxs_rnd(1:nFeat);
                    idxs_bbox_pair = idxAll(idxsPair,:);                    
                    feat = pt_get_temporal_features_img_dm(p, locs1, corres_pts1, locs2, corres_pts2, idxs_bbox_pair);
                    idxs = n_neg+1:n_neg+size(feat,1);
                    X_neg(idxs,:) = feat;
                    n_neg = n_neg + size(feat,1);
                end
            end
        end
    end
    
    check = 1;
end

% remove unused bins
X_pos(n_pos+1:end,:) = [];
X_neg(n_neg+1:end,:) = [];

mkdir_if_missing(p.ptPairwiseDir);
save(save_file, 'X_pos', 'X_neg', '-v7.3');

end
% ------------------------------------------------------------------------