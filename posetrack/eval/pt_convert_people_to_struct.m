function [ structPeople ] = pt_convert_people_to_struct( people, pidxs, parts )

num_joints = length(pidxs);

[numPersons, numFrames] = size(people);

structPeople = cell(numPersons,numFrames);

for p=1:numPersons
    for f=1:numFrames
        
        joints = people{p,f};
        if(isempty(joints))
            continue;
        end
        
        point = struct();
        
        np = 1;
        for j=1:num_joints
             pidx = pidxs(j);
             annopoint_idxs = parts(pidx+1).pos;
             assert(annopoint_idxs(1) == annopoint_idxs(2));
             
             if(any(isnan(joints(j,:))))
                 continue;
             end
             
             point(np).x = joints(j,1);
             point(np).y = joints(j,2);
             point(np).is_visible = true;
             point(np).id = annopoint_idxs(1);
             np = np + 1;
        end
        structPeople{p,f}.point = point;
    end
end


end

