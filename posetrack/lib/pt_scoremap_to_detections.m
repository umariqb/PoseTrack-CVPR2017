function [ detections ] = pt_scoremap_to_detections(p, score_fn, cidxs, scale)
%PT_SCOREMAP_TO_DETECTIONS Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 4)
    scale = 1;
end

stride              = p.stride;
half_stride         = stride/2;
locref_scale_mul    = p.locref_scale;
locref              = p.locref;
nextreg             = p.nextreg;
unLab_cls           = 'uint64';
max_sizet_val       = intmax(unLab_cls);

load(score_fn, 'scoremaps');
if locref
    load([score_fn, '_locreg'], 'locreg_pred');
end
if nextreg
    load([score_fn, '_nextreg'], 'nextreg_pred');
end

scales = p.scales;
num_scales = length(scales);
assert(num_scales == length(scoremaps)); 
locations = cell(num_scales,1);
I = cell(num_scales,1);
bVis = false;

if isfield(p, 'nms_dist')
    nms_dist = p.nms_dist;
else
    nms_dist = 1.5;
end

if(isfield(p, 'dets_per_part'))
    dets_per_part = p.dets_per_part;
else
    dets_per_part = 20;
end

detections = [];

for cidx = cidxs
    for s=1:num_scales
        p.scale_factor = scales(s)*scale;
        scale_factor = p.scale_factor;
        inv_scale_factor = 1/scale_factor;

        smaps = scoremaps{s};
        
        if(bVis)
            scmap_vis = visualise_scoremap(smaps);
            figure(119), imshow(scmap_vis);
            pause();
        end

        
        sm = smaps(:,:,cidx);
        locations = scoremap_to_detections(sm);
        I = nms_distance(locations, nms_dist);

        num_proposals = length(I);
        unPos = zeros(num_proposals, 2);
        unPos_sm = zeros(num_proposals, 2);
        if locref
            locationRefine = zeros(num_proposals, 1, 2);
        else
            locationRefine = [];
        end
        unProb = zeros(num_proposals, 1);
        if nextreg
            nextReg = zeros(num_proposals, size(nextreg_pred{1}, 3), size(nextreg_pred{1}, 4));
        else
            nextReg = [];
        end

        for k = 1:num_proposals
            idx = I(k);
            [row, col] = ind2sub(size(sm), idx);
            % transform heatmap to image coordinates
            %fprintf('row, col: %d %d\n', row, col);
            crd = [col-1, row-1]*stride;
            if p.res_net
               crd = crd + half_stride;
            end
            unPos(k, :)     = crd/scale_factor; %unnormalized position
            unPos_sm(k, :)  = [col, row]; % position in scoremaps
            unProb(k, :)    = sm(row,col); % probabilities
            if locref
               locreg_pred_s = locreg_pred{s};
               locationRefine(k, 1, :) = squeeze(locreg_pred_s(row, col, cidx, :))*locref_scale_mul*inv_scale_factor;
            end
            if nextreg
                nextreg_pred_s = nextreg_pred{s}*inv_scale_factor;
                nextReg(k, :, :) = squeeze(nextreg_pred_s(row, col, :, :));
            end
        end

        dets = MultiScaleDetections.make(unPos, unPos_sm, unProb, locationRefine, nextReg);
        dets.partClass = ones(size(unProb,1),1)*cidx;
        clearvars unPos unPos_sm unProb locationRefine nextReg;

        if p.nms_locref
            locRefJoint = squeeze(dets.locationRefine(:, 1, :));
            nd = size(dets.unProb,1);
            if(nd == 1)
                refinedCoord = dets.unPos + locRefJoint';
            else   
                refinedCoord = dets.unPos + locRefJoint;
            end
            locations = coord_to_scoremap(p, refinedCoord, [1 1]);

            if(bVis && false)
                imagesc(sm);
                figure(1);
                imagesc(sm);
                colorbar;

                % after refinement
                I = sub2ind(size(sm), locations(:,2)+1, locations(:,1)+1);
                mask = zeros(size(sm));
                mask(I) = 1;
                figure(2);
                imagesc(mask);

                % original
                I = sub2ind(size(sm), dets.unPos_sm(:,2), dets.unPos_sm(:,1));
                mask = zeros(size(sm));
                mask(I) = 1;
                figure(3);
                imagesc(mask);

                pause;


            end

            locations = [locations dets.unProb(:)];
            I = nms_distance(locations, p.nms_locref_dist);
%             fprintf('Scale %d: number of detections before %d after %d\n', s, size(dets.unProb, 1), length(I));
            dets = MultiScaleDetections.slice(dets, I);
        end
        dets.scale = ones(size(dets.unProb, 1), 1) * scales(s);
        if(isempty(detections))
            detections = dets;
        else
            detections = MultiScaleDetections.merge(detections, dets);
        end
    end
end

% if (isfield(p,'ignore_low_scores') && p.ignore_low_scores)
%     idxs = get_unary_idxs_local_sort_per_class(dets.unProb,dets_per_part);
%     fprintf('preserve %d/%d detections\n',length(idxs),length(dets.unProb));
%     detections = Detections.slice(detections, idxs);
% end

min_det_score = p.min_det_score;
detections.unProbNoThresh = detections.unProb;
detections.unProb(detections.unProb < min_det_score) = 0;

I = detections.unProb > min_det_score;
detections = MultiScaleDetections.slice(detections, I);

% perform nms accross scales
locRefJoint = squeeze(detections.locationRefine(:, 1, :));
refinedCoord = detections.unPos + locRefJoint;
patchSizes = p.patchSize * (1./detections.scale);


x1 = max(1,refinedCoord(:,1) - patchSizes/2);
x2 = refinedCoord(:,1) + patchSizes/2;
y1 = max(1,refinedCoord(:,2) - patchSizes/2);
y2 = refinedCoord(:,2) + patchSizes/2;

bbox = cat(2,x1,y1,x2,y2,detections.unProb);
I = nms_IoMin(bbox, 0.85);
detections = MultiScaleDetections.slice(detections, I);

