function [feat,idxs_bbox_pair,cb1,cb2] = pt_get_spatial_features_neighbour_img(locations,idxs_bbox_pair,next_joints, scales)
    
s = scales(idxs_bbox_pair(:,1));
s = [s,s];

cb1 = s.*locations(idxs_bbox_pair(:,1), :);
cb2 = s.*locations(idxs_bbox_pair(:,2), :);

deltaX = cb2(:,1)-cb1(:,1);
deltaY = cb2(:,2)-cb1(:,2);

nxt1 = squeeze(next_joints(idxs_bbox_pair(:,1), 1, :));
if(size(nxt1,2) == 1)
    nxt1 = nxt1';
end
deltaPred_forward = s.*nxt1;

nxt2 = squeeze(next_joints(idxs_bbox_pair(:,2), 2, :));
if(size(nxt2,2) == 1)
    nxt2 = nxt2';
end
deltaPred_backward = s.*nxt2;

if size(deltaPred_forward, 2) == 1
    deltaPred_forward = deltaPred_forward';
    deltaPred_backward = deltaPred_backward';
end

feat = cat(2, deltaX, deltaY, ...
              deltaPred_forward, deltaPred_backward);