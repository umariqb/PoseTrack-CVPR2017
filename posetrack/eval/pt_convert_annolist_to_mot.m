function annInfo = pt_convert_annolist_to_mot(annolist, isGT, pidxs, parts)


num_frames = annolist.num_frames;
num_persons = annolist.num_persons;
num_joints = length(pidxs);

annInfo.X = nan(num_frames, num_persons, num_joints);
annInfo.Y = nan(num_frames, num_persons, num_joints);
annInfo.isVis = nan(num_frames, num_persons, num_joints);

if(isGT)
    annInfo.rd    = nan(num_frames, num_persons, num_joints);
end

min_len = 1e10;
for f = 1:num_frames
    try
    for p=1:num_persons
        % take the info about which predictions should be considered for
        % evaluation
        if(isfield(annolist, 'keepPredInfo'))
            keepPredInfo = annolist.keepPredInfo{f};
        else
            keepPredInfo = ones(num_persons, num_joints);
        end
        [joints,visFlag] = pt_get_anno_joints_with_occ(annolist.annopoints{p,f}, pidxs, parts);
         for j=1:num_joints
            % if visible or set to zero due to occlusion flag
            if(visFlag(j) && keepPredInfo(p, j))
                annInfo.X(f,p,j) = joints(j,1);
                annInfo.Y(f,p,j) = joints(j,2);
                annInfo.isVis(f,p,j) = 1;
                if(isGT)
                    %reference distance
                    rect = annolist.annopoints{p, f}.head_rect;
                    BIAS = 1; %0.6 taken by MPII
                    annInfo.rd(f,p,j) = BIAS*norm(rect(1:2) - rect(3:4));
                    assert(~isnan(annInfo.rd(f,p,j)));
                end
            else
                annInfo.X(f,p,j) = 0;
                annInfo.Y(f,p,j) = 0;
                annInfo.isVis(f,p,j) = 0;
            end
        end
    end
    catch
        keyboard();
    end
end

end

