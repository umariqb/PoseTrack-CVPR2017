function [ joints, visFlags ] = pt_get_anno_joints_with_occ( annopoints, pidxs, parts )

num_joints = length(pidxs);
joints = NaN(num_joints, 2);
visFlags = zeros(num_joints,1);

if(~isfield(annopoints, 'point'))
    return;
end

points = annopoints.point;

if(~isfield(points, 'id'))
    return;
end

for j = 1:num_joints
    pidx = pidxs(j);
    annopoint_idxs = parts(pidx+1).pos;
    assert(annopoint_idxs(1) == annopoint_idxs(2));
    pt = util_get_annopoint_by_id(points, annopoint_idxs(1));
    if (~isempty(pt))
        joints(j, :) = [pt.x pt.y];
        if(isfield(pt, 'is_visible'))
            visFlags(j,1) = pt.is_visible;
        end
    end
end

end

