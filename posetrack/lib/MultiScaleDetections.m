classdef MultiScaleDetections
   
methods(Static)
    function res = make(unPos, unPos_sm, unProb, locationRefine, nextReg, unLab)
        res = struct();
        res.unPos = unPos;
        res.unPos_sm = unPos_sm;
        res.unProb = unProb;
        res.locationRefine = locationRefine;
        res.nextReg = nextReg;
        if nargin >= 6
            res.unLab = unLab;
        end
    end
    
    function dets = merge(dets1, dets2)
        dets = struct;
        dets.unPos = cat(1, dets1.unPos, dets2.unPos);
        dets.unProb = cat(1, dets1.unProb, dets2.unProb);
        if isfield(dets1, 'unPos_sm')
            dets.unPos_sm = cat(1, dets1.unPos_sm, dets2.unPos_sm);
        end
        if isfield(dets1, 'locationRefine')
            dets.locationRefine = cat(1, dets1.locationRefine, dets2.locationRefine);
        end
        if isfield(dets1, 'nextReg')
            dets.nextReg = cat(1, dets1.nextReg, dets2.nextReg);
        end
        if isfield(dets1, 'unProbNoThresh')
            dets.unProbNoThresh = cat(1, dets1.unProbNoThresh, dets2.unProbNoThresh);
        end
        if isfield(dets1, 'unLab')
            dets.unLab = cat(1, dets1.unLab, dets2.unLab);
        end
        if isfield(dets1, 'scale')
            dets.scale = cat(1, dets1.scale, dets2.scale);
        end
        if isfield(dets1, 'frameIndex')
            dets.frameIndex = cat(1, dets1.frameIndex, dets2.frameIndex);
        end
        if isfield(dets1, 'index')
            dets.index = cat(1, dets1.index, dets2.index);
        end
        if isfield(dets1, 'partClass')
            dets.partClass = cat(1, dets1.partClass, dets2.partClass);
        end


    end

    function dets = slice(dets, idxs)
        dets.unPos = dets.unPos(idxs, :);
        dets.unProb = dets.unProb(idxs, :);

        if isfield(dets, 'unPos_sm')
            if ~isempty(dets.unPos_sm)
                dets.unPos_sm = dets.unPos_sm(idxs, :);
            end
        end
        
        if isfield(dets, 'locationRefine')
            if ~isempty(dets.locationRefine)
                dets.locationRefine = dets.locationRefine(idxs, :, :);
            end
        end
        if isfield(dets, 'nextReg')
            if ~isempty(dets.nextReg)
                dets.nextReg = dets.nextReg(idxs, :, :);
            end
        end
        if isfield(dets, 'unProbNoThresh')
            dets.unProbNoThresh = dets.unProbNoThresh(idxs, :);
        end
        if isfield(dets, 'unLab') && ~isempty(dets.unLab)
            dets.unLab = dets.unLab(idxs, :);
        end
        if isfield(dets, 'scale') && ~isempty(dets.scale)
            dets.scale = dets.scale(idxs, :);
        end
        if isfield(dets, 'frameIndex') && ~isempty(dets.frameIndex)
            dets.frameIndex = dets.frameIndex(idxs, :);
        end
        if isfield(dets, 'index') && ~isempty(dets.index)
            dets.index = dets.index(idxs, :);
        end
        if isfield(dets, 'partClass') && ~isempty(dets.partClass)
            dets.partClass = dets.partClass(idxs, :);
        end
        if isfield(dets, 'Pose') && ~isempty(dets.Pose)
            dets.Pose = dets.Pose(idxs, :);
        end

    end
    
    function save(out_dets, unLab, fname)
        unPos = out_dets.unPos;
        unProb = out_dets.unProb;
        if isfield(out_dets, 'locationRefine')
            locationRefine = out_dets.locationRefine;
        end
        %unProb = unProbNoThresh;
        save(fname,'unLab','unPos','unProb','locationRefine');
    end
end

end