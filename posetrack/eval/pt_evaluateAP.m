function [apAll] = pt_evaluateAP(expidxs, options)
% implementation of AP measure,
% as defined in [Pishchulin et al., arXiv'15]

if(nargin < 2)
    options.td = 0.2;
    options.minTrackLen = 5;
    options.minAvgJoints = 3;
end


partNames = {'right ankle','right knee','right hip','left hip','left knee','left ankle','right wrist','right elbow','right shoulder','left shoulder','left elbow','left wrist','neck','top head','avg full body'};
colorIdxs = [1 1];

fprintf('pt_evaluateAP()\n');

thresh = options.td;

tableTex = cell(length(expidxs)+1,1);

% load ground truth
p = pt_exp_params(1);
gtAnnolist = load(p.testGT,'annolist');
gtAnnolist = gtAnnolist.annolist;
parts = get_parts();
pidxs = p.pidxs;
gtAnnolistMPII = pt_convert_to_mpii_format(gtAnnolist);

apAll = zeros(15,length(expidxs));
markerAll = repmat({'-'},length(expidxs),1);
% tableDir = [fileparts(p.predFilename) '/latex']; if (~exist(tableDir,'dir')), mkdir(tableDir); end

for i = 1:length(expidxs);
    
    expidx = expidxs(i);
    p = pt_exp_params(expidx);
    pt_generate_result_file(expidx, options);
    
    num_joints = length(pidxs);

    annPath = fullfile(p.ptResultsDir, ['exp', num2str(expidx)], 'pred_annolist');
    estAnnolist = load(annPath, 'annolist');
    estAnnolist = estAnnolist.annolist;
    estAnnolistMPII = pt_convert_to_mpii_format(estAnnolist);
    
    assert(length(gtAnnolistMPII) == length(estAnnolistMPII));
    % assign predicted poses to GT poses
    [scoresAll, labelsAll, nGTall] = pt_assignGTmulti(estAnnolistMPII,gtAnnolistMPII,thresh(end));
    % compute average precision (AP) per part
    ap = zeros(size(nGTall,1)+1,1);
    for j = 1:size(nGTall,1)
      scores = []; labels = [];
      for imgidx = 1:length(gtAnnolistMPII)
        scores = [scores; scoresAll{j}{imgidx}];
        labels = [labels; labelsAll{j}{imgidx}];
      end
      % compute precision/recall
      [precision,recall] = getRPC(scores,labels,sum(nGTall(j,:)));
      % compute AP
      ap(j) = VOCap(recall,precision)*100;
    end
    % compute mean AP
    ap(end) = mean(ap(1:end-1));
    columnNames = partNames;
%     save([fileparts(p.predFilename) '/apAll'],'ap','columnNames');
    % plot results
    [row, header] = pt_genTableAP(ap,p.name);
    tableTex{1} = header;
    tableTex{i+1} = row;
    apAll(:,i) = ap;
end

% fid = fopen([tableDir '/ap.tex'],'wt');assert(fid ~= -1);
% for i=1:length(tableTex),fprintf(fid,'%s',tableTex{i}); end; fclose(fid);
% 
end