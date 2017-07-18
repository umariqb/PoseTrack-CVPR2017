function p = pt_exp_params(expidx)

p = [];
p.code_dir = pwd();
p.expDir = [p.code_dir '/data/'];
p.vidDir = [p.code_dir '/data/videos/'];
p.latexDir = [p.code_dir '/output/latex/'];
p.plotsDir = [p.code_dir '/output/plots/'];
p.ptSolver  = [p.code_dir '/posetrack/solver/pt-solver-callback'];
p.gurobi_license_file = [p.expDir '/gurobi/gurobi.lic;'];

p.deepMatching = [p.code_dir '/external/deepmatching/deepmatching'];
p.pairwise_relations = [p.code_dir '/deepcut/data/pairwise/all_pairs_stats_all.mat'];
p.liblinear_predict = true;

switch expidx
    case 1
        p.name = 'Bonn multiperson posetrack test';
        p.shortName = 'bonn-multiperson-posetrack';
        p.vidDir = fullfile(p.expDir, p.shortName, 'videos');
        
        p.testGT = fullfile(p.expDir, p.shortName, 'annolist/test/annolist');
        p.testPad = 0;

        p.trainGT = fullfile(p.expDir, p.shortName, 'annolist/train/annolist');
        p.trainPad = 0;
        
        p.evalTest = fullfile(p.expDir, p.shortName, 'test', 'predictions');
        p.outputDir = [p.expDir '/' p.shortName '/output/'];
        p.ptMulticutDir = [p.expDir '/' p.shortName '/pt-multicut/'];
        p.ptDetectionsDir = [p.expDir '/' p.shortName '/detections/'];
        p.ptResultsDir = [p.expDir '/' p.shortName '/results/'];

        p.pairwiseDir = fullfile(p.expDir, 'pairwise');
        p.ptPairwiseDir = fullfile(p.expDir, 'pt-pairwise');
        
        p.net_dir = fullfile(p.code_dir, 'deepcut/data/caffe-models');
        p.net_def_file = 'ResNet-101-FCN_out_14_sigmoid_locreg_allpairs_test.prototxt';
        p.net_bin_file = 'ResNet-101-mpii-multiperson.caffemodel';

        p.unary_scoremap_dir = fullfile(p.expDir, p.shortName, 'scoremaps', 'test');
        p.pairwise_scoremap_dir = fullfile(p.expDir, p.shortName, 'scoremaps', 'test');
        p.scoremaps = fullfile(p.expDir, p.shortName, 'scoremaps');
        p.correspondences = fullfile(p.expDir, p.shortName, 'correspondences');
        p.corresName = '%05d_%05d.txt';
        p.stride = 8;
        p.scale_factor = 1;
        p.multiscale = true;
        p.scales = 0.6:0.3:1.5;
        p.locref = true;
        p.nextreg = true;
        p.allpairs = true;
        p.patchSize   = 70;

	    p.res_net = true;

        p.nms_dist = 7.5;
        p.nms_locref = true;
        p.nms_locref_dist = 7.5;

        p.pidxs = [0 2 4 5 7 9 12 14 16 17 19 21 22 23];
        p.cidxs_full = 1:14;
        p.temporal_cidxs = [13];

        p.stagewise = true;
        p.split_threshold = 0.4;
        p.cidxs_stages = {[9 10 13 14], [9 10 13 14 7 8 11 12], [9 10 13 14 7 8 11 12 1 2 3 4 5 6]};
        p.correction_stages = [0 0 0];

        p.all_parts_on = false;
        p.nFeatSample = 10^2;

        p.dets_per_part = 20;

        p.ignore_low_scores = true;
        p.min_det_score = 0.2;
        p.min_det_score_face = 0.2;
        p.single_people_solver = false;
        p.high_scores_per_class = true;
        p.all_parts_on = false;
        p.multi_people = true;
        p.time_limit = 86400;
        
        p.scale = 4;
        p.colorIdxs = [5 1];
        p.refHeight = 400;

        p.multicut = true;
        
        p.maxFrameDist = 3; 
        p.maxDim   = 2000;
        p.maxDimDM = 1440;
        p.matchRadius = 100;
        
        p.temporalWinSize = 31;
        p.trackMinLen = 0;
        
    case 2
        p = pt_exp_params(1);
        p.name = 'Exp2: Head Only';
        p.temporal_cidxs = [13]; % head only
    case 3
        p = pt_exp_params(1);
        p.name = 'Exp3:  head, neck, shoulders, Frames = 3\n';
        p.temporal_cidxs = [14,13,9,10]; % head, neck, shoulders
    case 4
        p = pt_exp_params(1);
        p.name = 'Exp4: head, neck, shoulders, hips';
        p.temporal_cidxs = [14,13,9,10,3,4]; % head, neck, shoulders, hips
    case 5
        p = pt_exp_params(1);
        p.name = 'Exp5: head, wrists, ankles';
        p.temporal_cidxs = [14,7,12,1,6]; % head, wrists, ankles
    case 6
        p = pt_exp_params(1);
        p.name = 'Exp6: Head Only, Frames = 5';
        p.temporal_cidxs = [13]; % head only
        p.maxFrameDist = 5;
    case 7
        p = pt_exp_params(1);
        p.name = 'Exp7:  head, neck, shoulders. Frames = 5';
        p.temporal_cidxs = [14,13,9,10]; % head, neck, shoulders
        p.maxFrameDist = 5;
    case 8
        p = pt_exp_params(1);
        p.name = 'Exp8: Head Only, Frames = 7';
        p.temporal_cidxs = [13]; % head only
        p.maxFrameDist = 7;
    case 9
        p = pt_exp_params(1);
        p.name = 'Exp9:  head, neck, shoulders. Frames = 7';
        p.temporal_cidxs = [14,13,9,10]; % head, neck, shoulders
        p.maxFrameDist = 7;
    case 10
        p = pt_exp_params(1);
        p.name = 'Exp10: Head Only. Frames = 1';
        p.temporal_cidxs = [13]; % head only
        p.maxFrameDist = 1;
    case 11
        p = pt_exp_params(1);
        p.name = 'Exp11:  head, neck, shoulders. Frames = 1\n';
        p.temporal_cidxs = [14,13,9,10]; % head, neck, shoulders
        p.maxFrameDist = 1;
    case 12
        p = pt_exp_params(1);
        p.name = 'Exp12: head, neck, shoulders, hips. Frames = 1';
        p.temporal_cidxs = [14,13,9,10,3,4]; % head, neck, shoulders, hips
        p.maxFrameDist = 1;
    case 13
        p = pt_exp_params(1);
        p.name = 'Exp13: head, wrists, ankles. Frames = 1';
        p.temporal_cidxs = [14,7,12,1,6]; % head, wrists, ankles
        p.maxFrameDist = 1;
end

if (isfield(p,'colorIdxs') && ~isempty(p.colorIdxs))
    p.colorName = eval_get_color_new(p.colorIdxs);
    p.colorName = p.colorName ./ 255;
end

p.exp_dir = fullfile(p.expDir, p.shortName);
p.res_net = isfield(p, 'res_net') && p.res_net;
p.rpn_detect = isfield(p, 'rpn_detect') && p.rpn_detect;
p.detcrop_recall = isfield(p, 'detcrop_recall') && p.detcrop_recall;
p.detcrop_image = isfield(p, 'detcrop_image') && p.detcrop_image;
p.histogram_pairwise = isfield(p, 'histogram_pairwise') && p.histogram_pairwise;
p.nms_locref = isfield(p, 'nms_locref') && p.nms_locref;
p.stagewise = isfield(p, 'stagewise') && p.stagewise;
if p.stagewise
    p.num_stages = length(p.cidxs_stages);
end
p.person_part = isfield(p, 'person_part') && p.person_part;
p.complete_clusters = isfield(p, 'complete_clusters') && p.complete_clusters;

p.locref_scale = sqrt(53);

p.mean_pixel = [104, 117, 123];

if ~isfield(p, 'stride')
    p.stride = 8;
end



end
