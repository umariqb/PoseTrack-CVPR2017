function [expMetrics] = pt_eval_pose_tracking( expidxs, options, bRecompute)

fprintf('pt_eval_pose_tracking()\n');

if(nargin < 2)
    fprintf('Using default options.\n');
    options.td = 0.2;
    options.minTrackLen = 5;
    options.minAvgJoints = 3;
end

if(nargin < 3)
    bRecompute = true;
end


expMetrics = cell(size(expidxs));
p = pt_exp_params(1);
gtAnnolist = load(p.testGT,'annolist');
gtAnnolist = gtAnnolist.annolist;
parts = get_parts();
pidxs = p.pidxs;

headerPrinted = false;
for expidx = expidxs
    p = pt_exp_params(expidx);
    if(bRecompute)
        pt_generate_result_file(expidx, options);
    end
    
    num_joints = length(pidxs);

    annPath = fullfile(p.ptResultsDir, ['exp', num2str(expidx)], 'pred_annolist');
    estAnnolist = load(annPath, 'annolist');
    estAnnolist = estAnnolist.annolist;
    
    assert(length(estAnnolist) == length(gtAnnolist));
    num_videos = length(gtAnnolist);
    
    metricsSum = zeros(num_joints,12); 
    for vIdx = 1:num_videos
        gtAnn = gtAnnolist(vIdx);
        estAnn = estAnnolist(vIdx);
        
        gtAnnMPII  = pt_convert_to_mpii_format(gtAnn);
        estAnnMPII = pt_convert_to_mpii_format(estAnn);
        assert(length(gtAnnMPII) == length(estAnnMPII));

        % get the information of prediction that should not be considered
        % during evaluation of MOT due to occlusion flags
        [~, ~, ~, keepPredInfo] = pt_assignGTmulti(estAnnMPII,gtAnnMPII,options.td);
        estAnn = pt_assign_keep_pred_info(estAnn, keepPredInfo);
        assert(strcmp(gtAnn.name, estAnn.name));
        gtInfo = pt_convert_annolist_to_mot(gtAnn,true,pidxs,parts);
        assert(isfield(gtInfo, 'isVis') && isfield(gtInfo, 'rd'));
        stateInfo = pt_convert_annolist_to_mot(estAnn,false,pidxs,parts);
        
        for j = 1:num_joints
            gtInfoJ.X = squeeze(gtInfo.X(:,:,j));
            gtInfoJ.Y = squeeze(gtInfo.Y(:,:,j));
            gtInfoJ.isVis = squeeze(gtInfo.isVis(:,:,j));
            gtInfoJ.rd = squeeze(gtInfo.rd(:,:,j));

            stateInfoJ.X = squeeze(stateInfo.X(:,:,j));
            stateInfoJ.Y = squeeze(stateInfo.Y(:,:,j));
            stateInfoJ.isVis = squeeze(stateInfo.isVis(:,:,j));
            [m, info, ~] = pt_CLEAR_MOT_HUN(gtInfoJ, stateInfoJ, options);
            % sum all metrics for each joint
            metricsSum(j,:) = metricsSum(j,:) + m;
        end
    end

    % compute MOTA and AP for each joint
    % Order: [1-Ngt, 2-MT, 3-PT, 4-ML, 5-falsepositives, 6-missed, 7-idswitches, 8-FRA, 9-numGT, 10-numC, 11-numD, 12-Fgt]
    falsepositives  = metricsSum(:,5);
    missed          = metricsSum(:,6);
    idswitches      = metricsSum(:,7);
    numGT           = metricsSum(:,9);
    numC            = metricsSum(:,10);

    MOTA = (ones(num_joints,1) - ((falsepositives + missed + idswitches)./ numGT)) * 100;
    precision   = numC ./ (falsepositives + numC) * 100;

    perJointMOTA = pt_compute_joint_metrics_to_print(MOTA);
    perJointPrecision = pt_compute_joint_metrics_to_print(precision);
    perJointMetrics = [perJointPrecision;perJointMOTA];

    % compute all metrics for complete dataset (including all joints)
    metricsSumDS    = sum(metricsSum);
    mc = num2cell(metricsSumDS);
    [Ngt, MT, PT, ML, falsepositives, missed, idswitches, FRA, numGT, numC, numD, Fgt] = deal(mc{:});
    td = options.td;
    MOTP=(1-numD/numC/td) * 100; % avg distance to [0,100] 
    MOTA = (1 - ((falsepositives + missed + idswitches)/ numGT)) * 100;
    recall=numC/numGT*100;
    precision = numC ./ (falsepositives + numC) * 100;
    FAR=falsepositives/Fgt;
    
    metrics = [recall, precision, FAR, Ngt, MT, PT, ML, falsepositives, ...
        missed, idswitches, FRA, MOTA, MOTP];
        
    dispMetrics = [1,2,4,5,6,7,10,11,12,13];
    if(~headerPrinted)
        pt_printMetrics(metrics, perJointMetrics, false, dispMetrics);
        headerPrinted = true;
    else
        pt_printMetrics(metrics, perJointMetrics, false, dispMetrics);
    end
    expMetrics{expidx} = metrics;
        
end

end

function outMetrics = pt_compute_joint_metrics_to_print(inMetrics)

head = [13,14]; 
sho  = [9,10];  
elb  = [8,11];  
wri  = [7,12];  
hips = [3,14];
knee = [2,5];
ank  = [1,6];

outMetrics = zeros(1,7);
outMetrics(1) = mean(inMetrics(head));
outMetrics(2) = mean(inMetrics(sho));
outMetrics(3) = mean(inMetrics(elb));
outMetrics(4) = mean(inMetrics(wri));
outMetrics(5) = mean(inMetrics(hips));
outMetrics(6) = mean(inMetrics(knee));
outMetrics(7) = mean(inMetrics(ank));

end

function [ann] = pt_assign_keep_pred_info(ann, keepPredInfo)
    for f=1:ann.num_frames
        keepInfo = [];
        idxs = ~cellfun(@isempty,ann.annopoints(:,f));
        keepInfo(idxs,:) =  keepPredInfo{f};
        ann.keepPredInfo{f} = keepInfo;
    end
end
        

