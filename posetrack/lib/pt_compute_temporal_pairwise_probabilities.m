function [ pwProb ] = pt_compute_temporal_pairwise_probabilities(p, detections, temporal_model, cidx, fn, corres_dir)
%PT_COMPUTE_PAIRWISE_PROBABILITIES Summary of this function goes here
%   Detailed explanation goes here

if (ischar(cidx))
    cidx = str2num(cidx);
end

pw_weight = 1;

frames = min(detections.frameIndex):max(detections.frameIndex);
num_frames = length(frames);

nSpatDetRel = 0;
spatial_det_rel = cell(num_frames,1);
for f = 1:num_frames
    fIdx = frames(f);
    num_dets  = sum(detections.frameIndex == fIdx);
    spatial_det_rel{f} = pt_build_frame_pairs(num_dets, num_dets); %FIXME: change the name of the function to make it general.
    nSpatDetRel  = nSpatDetRel + size(spatial_det_rel{f}, 1); 
end

frame_pairs = pt_build_frame_pairs(num_frames, p.maxFrameDist);
temporal_det_rel = cell(size(frame_pairs,1),1);
nTempDetRel = 0;
for f = 1:length(temporal_det_rel)
    pair        = frames(frame_pairs(f,:));
    dets1 = MultiScaleDetections.slice(detections, detections.frameIndex == pair(1));
    dets2 = MultiScaleDetections.slice(detections, detections.frameIndex == pair(2));
    
    num_dets1   = size(dets1.unProb, 1);
    num_dets2   = size(dets2.unProb, 1);
    [q1,q2]     = meshgrid(1:num_dets1,1:num_dets2);
    idxsAllrel  = [q1(:) q2(:)];
    temporal_det_rel{f} = idxsAllrel;
    nTempDetRel = nTempDetRel + size(idxsAllrel, 1);
end

pwProb = nan(nSpatDetRel+nTempDetRel, 6);

idxStart = 1;

% probabilities for with-in frame detections
for f = 1:num_frames
    fIdx = frames(f);
    dets = MultiScaleDetections.slice(detections, detections.frameIndex == fIdx);
    dRel  = spatial_det_rel{f};
    num_dets = size(dets.unProb,1);
    
    if(num_dets == 0)
        continue;
    end
    
    if p.locref
        try
        if(num_dets == 1)
            locations = dets.unPos+(squeeze(dets.locationRefine(:,1,:))');
        else   
            locations = dets.unPos+squeeze(dets.locationRefine(:,1,:));
        end
        catch
            keyboard();
        end
    else
        locations = dets.unPos;
    end
    
    if(size(locations,2) == 2)
        scales = dets.scale;
        ps = p.patchSize * 1./scales;
        r  = ps/2;
        x1 = locations(:,1) - r;
        y1 = locations(:,2) - r ;
        x2 = x1 + ps;
        y2 = y1 + ps; 
        bbox = cat(2, x1,y1,x2,y2);
    else
        bbox = locations;
    end
    
    IoU = boxoverlap_one2one(bbox(dRel(:,1),:), bbox(dRel(:,2),:));
    
    if(sum(isnan(IoU)))
        keyboard();
    end
    
    idxs = idxStart:idxStart+size(dRel,1)-1;
    pwProb(idxs,1:2) = dets.frameIndex(1);
    pwProb(idxs,3) = dets.index(dRel(:,1));
    pwProb(idxs,4) = dets.index(dRel(:,2));
    pwProb(idxs,5) = cidx;
    pwProb(idxs,6) = IoU(:,1);

    idxStart = idxStart + size(dRel,1);
end

% probabilities for accross frame detections
for f = 1:size(frame_pairs,1)    
    
    fPair = frames(frame_pairs(f,:));
    dRel  = temporal_det_rel{f};

    fr_fn1 = fn(fPair(1)).name;
    fr_fn2 = fn(fPair(2)).name;

    [~,fr_name1,~] = fileparts(fr_fn1);
    [~,fr_name2,~] = fileparts(fr_fn2);

    corres_fn = fullfile(corres_dir, [fr_name1,'_',fr_name2,'.txt']);
    corres_pts = pt_load_dm_correspondences(corres_fn);
    corres_pts1 = corres_pts(:,1:2);
    corres_pts2 = corres_pts(:,3:4);
   
    dets1 = MultiScaleDetections.slice(detections, detections.frameIndex == fPair(1));
    dets2 = MultiScaleDetections.slice(detections, detections.frameIndex == fPair(2));
    
    num_dets1 = size(dets1.unProb,1);
    num_dets2 = size(dets2.unProb,1);
    
    if(num_dets1 == 0 || num_dets2 == 0)
        continue;
    end

    locs1 = cat(2, dets1.unPos, dets1.scale, dets1.unProb);
    locs2 = cat(2, dets2.unPos, dets2.scale, dets2.unProb);
    
    feat = pt_get_temporal_features_img_dm(p, locs1, corres_pts1, locs2, corres_pts2, dRel);
%     feat = feat(:,1:5);
    feat_norm = getFeatNorm(feat,temporal_model.diff.training_opts.X_min,temporal_model.diff.training_opts.X_max);
    ex = sparse(double(feat_norm));
    model = temporal_model.diff.log_reg;

    if p.liblinear_predict
        [~,~,prob_out] = predict(zeros(size(ex,1),1), ex, model, '-b 1 -q');
        prob = prob_out(:,1);
    else
        x = [ full(ex) ones(size(ex, 1), 1)];
        prob = 1./(1+exp(-x*model.w'));
    end
        
    prob = prob.^pw_weight;
    
    if(sum(isnan(prob)))
        keyboard();
    end

    idxs = idxStart:idxStart+size(dRel,1)-1;
    pwProb(idxs,1) = fPair(1);
    pwProb(idxs,2) = fPair(2);
    pwProb(idxs,3) = dets1.index(dRel(:,1));
    pwProb(idxs,4) = dets2.index(dRel(:,2));
    pwProb(idxs,5) = cidx;
    pwProb(idxs,6) = prob(:,1);

    idxStart = idxStart + size(dRel,1);
end

pwProb(idxStart:end,:) = [];

end

