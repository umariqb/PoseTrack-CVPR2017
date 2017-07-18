function [annolist] = pt_convert_to_mpii_format(inAnnolist)

total_frames = sum([inAnnolist.num_frames]);
annolist(1:total_frames) = struct();
num_videos = length(inAnnolist);
nimgs = 1;
for vIdx = 1:num_videos
    num_frames = inAnnolist(vIdx).num_frames;
    num_persons = inAnnolist(vIdx).num_persons;
    
    for fIdx=1:num_frames
       npersons = 1;
       for pIdx=1:num_persons
            ann = inAnnolist(vIdx).annopoints{pIdx,fIdx};
            if(isempty(ann))
                continue;
            else
                annolist(nimgs).annorect(npersons).annopoints.point = ann.point;
                if(isfield(ann, 'head_rect'))
                    annolist(nimgs).annorect(npersons).x1 = ann.head_rect(1);
                    annolist(nimgs).annorect(npersons).y1 = ann.head_rect(2);
                    annolist(nimgs).annorect(npersons).x2 = ann.head_rect(3);
                    annolist(nimgs).annorect(npersons).y2 = ann.head_rect(4);
                end
            end
            npersons = npersons + 1;
       end
       nimgs=nimgs+1;
    end
end

end