function pt_generate_result_file( expidx, options)

if(nargin < 2)
    fprintf('Using default options.\n');
    options.minTrackLen = 5;
    options.minAvgJoints = 3;
end

p = pt_exp_params(expidx);
load(p.testGT,'annolist');
gtAnnolist = annolist;
parts = get_parts();

num_videos = length(gtAnnolist);

annolist = struct();


for vIdx = 1:num_videos
    
    vid_name = gtAnnolist(vIdx).name;
    gtNumFrames = gtAnnolist(vIdx).num_frames;
%     p.ptMulticutDir = '/home/ibal_109/work/2016/posetrack/data/bonn-multiperson-posetrack/pt-multicut_final';
    predPath = [p.ptMulticutDir '/prediction_' num2str(expidx) '_' vid_name, '.mat'];
    if(exist(predPath, 'file'))
        load(predPath, 'people');
        
        numPeople = size(people,1);
        numFrames = size(people,2);
        
%         assert(numFrames == gtNumFrames, sprintf('Number of frames not equal for video %d. (%d,%d)\n', vIdx, numFrames, gtNumFrames));
        corrIdxs = zeros(1,numPeople);
        for pIdx=1:size(people,1)
            len = sum(~cellfun(@isempty,people(pIdx,:)));
            if(len > options.minTrackLen)%p.trackMinLen)
                corrIdxs(pIdx) = 1;
            end
        end
        people = people(logical(corrIdxs), :);
        numPeople = size(people,1);
        corrIdxs = zeros(1,numPeople);
        for pIdx=1:size(people,1)
            nj = 0;
            nf = 0;
            for fIdx=1:size(people,2)
                joints = people{pIdx,fIdx};
                if(isempty(joints))
                    continue;
                end
                num = sum(~isnan(joints(:)))/2;
                nj = nj + num;
                nf = nf + 1;
            end
            avg = nj / nf;
            if(avg > options.minAvgJoints)
                corrIdxs(pIdx) = 1;
            end
        end

        people = people(logical(corrIdxs), :);
        
        annopoints = pt_convert_people_to_struct(people, p.pidxs, parts);
        annolist(vIdx).num_frames = size(annopoints, 2);
        annolist(vIdx).name = vid_name;
        annolist(vIdx).num_persons = size(annopoints, 1);
        annolist(vIdx).annopoints = annopoints;
    else
        fprintf('Video %d: Prediction does not exist\n', vIdx);
    end
end

savePath = [p.ptResultsDir, 'exp', num2str(expidx)];
mkdir_if_missing(savePath);
fprintf('Saving annolist file at: %s\n', savePath);
save([savePath,'/pred_annolist'], 'annolist');

end

