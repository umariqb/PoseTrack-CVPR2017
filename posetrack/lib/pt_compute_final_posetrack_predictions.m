function [ people ] = pt_compute_final_posetrack_predictions( dets )

num_joints = 14;

unLab = dets.unLab;

people = cell(1, 1);

frames = [min(dets.frameIndex):max(dets.frameIndex)];

for fIdx = frames
    
    fidxs = dets.frameIndex == fIdx;
    f_unLab = unLab(fidxs,:);
    fDetections = MultiScaleDetections.slice(dets, fidxs);
    flc_uniq = unique(f_unLab(:,2));

    for j = 1:length(flc_uniq)
        lc = flc_uniq(j);
        idxsc = find(f_unLab(:,2) == lc);
        lp_uniq = unique(fDetections.partClass(idxsc));

        keypoints = nan(num_joints,2);
        for i = 1:length(lp_uniq)
            lp = lp_uniq(i);
            if lp > 1000
                continue;
            end
            idxsp = find(fDetections.partClass(idxsc,1) == lp);

            prob = dets.unProb(idxsc(idxsp));
            w = prob;
            w = w./sum(w);

            loc_refine = squeeze(fDetections.locationRefine(idxsc(idxsp), 1, :));
            if size(loc_refine, 2) == 1
                loc_refine = loc_refine';
            end
            pos = sum((fDetections.unPos(idxsc(idxsp),:) + loc_refine).*[w w],1);

            keypoints(lp, :) = pos;
        end
    
        people{lc+1,fIdx} = keypoints;
    end
end

