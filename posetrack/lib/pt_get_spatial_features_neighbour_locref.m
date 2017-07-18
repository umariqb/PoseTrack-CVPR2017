function [feat,idxs_bbox_pair,cb1,cb2] = pt_get_spatial_features_neighbour_locref(locations,idxs_bbox_pair,next_joints, locref, scales)

s = scales(idxs_bbox_pair(:,1));
s = [s,s];

cb1 = s.*locations(idxs_bbox_pair(:,1), :);
cb2 = s.*locations(idxs_bbox_pair(:,2), :);

delta = cb2 - cb1;

deltaReal_forward = delta + squeeze(s.*locref(idxs_bbox_pair(:,2), 2, :));
deltaReal_backward = -delta + squeeze(s.*locref(idxs_bbox_pair(:,1), 1, :));

deltaPred_forward = squeeze(s.*next_joints(idxs_bbox_pair(:,1), 1, :));
deltaPred_backward = squeeze(s.*next_joints(idxs_bbox_pair(:,2), 2, :));

feat = cat(2, deltaReal_forward, deltaReal_backward, ...
              deltaPred_forward, deltaPred_backward);
