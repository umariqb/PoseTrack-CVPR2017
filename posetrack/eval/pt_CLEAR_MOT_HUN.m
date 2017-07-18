function [metrics metricsInfo additionalInfo]=CLEAR_MOT_HUN(gtInfo,stateInfo, options)
% compute CLEAR MOT and other metrics
%
% metrics contains the following
% [1]   recall	- recall = percentage of detected targets
% [2]   precision	- precision = percentage of correctly detected targets
% [3]   FAR		- number of false alarms per frame
% [4]   GT        - number of ground truth trajectories
% [5-7] MT, PT, ML	- number of mostly tracked, partially tracked and mostly lost trajectories
% [8]   falsepositives- number of false positives (FP)
% [9]   missed        - number of missed targets (FN)
% [10]  idswitches	- number of id switches     (IDs)
% [11]  FRA       - number of fragmentations
% [12]  MOTA	- Multi-object tracking accuracy in [0,100]
% [13]  MOTP	- Multi-object tracking precision in [0,100] (3D) / [td,100] (2D)
% [14]  MOTAL	- Multi-object tracking accuracy in [0,100] with log10(idswitches)
%
% 
% (C) Anton Milan, 2012-2014

if(nargin < 3)
    options = struct();
end

if ~isfield(options,'td')
    options.td=0.2;
end

[Fgt, Ngt]=size(gtInfo.X);
[F, N]=size(stateInfo.X);
td=options.td;

% if stateInfo shorter, pad with zeros
if F<Fgt
    missingFrames = F+1:Fgt;
    stateInfo.X(missingFrames,:,:)=0;
    stateInfo.Y(missingFrames,:,:)=0;
    stateInfo.isVis(missingFrames,:,:)=0;
end
[F, N]=size(stateInfo.X);

%  remove info of the targets that are always occluded
visGTInd = ~~sum(gtInfo.X,1);
gtInfo.X = gtInfo.X(:,visGTInd);
gtInfo.Y = gtInfo.Y(:,visGTInd);
gtInfo.rd = gtInfo.rd(:,visGTInd);
gtInfo.isVis=gtInfo.isVis(:,visGTInd);
[Fgt, Ngt]=size(gtInfo.X);
gtInd=~~gtInfo.X;

visStateInd = ~~sum(stateInfo.X,1);
stateInfo.X = stateInfo.X(:,visStateInd);
stateInfo.Y = stateInfo.Y(:,visStateInd);
stateInfo.isVis=stateInfo.isVis(:,visStateInd);
[F, N]=size(stateInfo.X);
stInd=~~stateInfo.X;

%%
metricsInfo.names.long = {'Recall','Precision','False Alarm Rate', ...
    'GT Tracks','Mostly Tracked','Partially Tracked','Mostly Lost', ...
    'False Positives', 'False Negatives', 'ID Switches', 'Fragmentations', ...
    'MOTA','MOTP', 'MOTA Log'};

metricsInfo.names.short = {'Rcll','Prcn','FAR', ...
    'GT','MT','PT','ML', ...
    'FP', 'FN', 'IDs', 'FM', ...
    'MOTA','MOTP', 'MOTAL'};

metricsInfo.widths.long = [6 9 16 9 14 17 11 15 15 11 14 5 5 8];
metricsInfo.widths.short = [5 5 5 3 3 3 3 4 4 3 3 5 5 5];

metricsInfo.format.long = {'.1f','.1f','.2f', ...
    'i','i','i','i', ...
    'i','i','i','i', ...
    '.1f','.1f','.1f'};

metricsInfo.format.short=metricsInfo.format.long;


metrics=zeros(1,12);
metrics(9)=numel(find(gtInd));  % False Negatives (missed)
metrics(7)=Ngt;                 % Mostly Lost
metrics(4)=Ngt;                 % GT Trajectories

additionalInfo=[];
% nothing to be done, if state is empty
if ~N, return; end

% mapping
M=zeros(F,Ngt);

mme=zeros(1,F); % ID Switchtes (mismatches)
c=zeros(1,F);   % matches found
fp=zeros(1,F);  % false positives
m=zeros(1,F);   % misses = false negatives
g=zeros(1,F);
d=zeros(F,Ngt);  % all distances;
ious=Inf*ones(F,Ngt);  % all overlaps

alltracked=zeros(F,Ngt);
allfalsepos=zeros(F,N);

g = [];
for t=1:Fgt
    g(t)=numel(find(gtInd(t,:)));
    
    % mapping for current frame
    if t>1   % 
        mappings=find(M(t-1,:));
        for map=mappings                              % finds matches between first frame and second frame
            if gtInd(t,map) && stInd(t,M(t-1,map)) && matched(gtInfo,stateInfo,t,map,M(t-1,map),td)
                M(t,map)=M(t-1,map);
            end
        end
    end
    
    GTsNotMapped=find(~M(t,:) & gtInd(t,:)); 
    EsNotMapped=setdiff(find(stInd(t,:)),M(t,:));

    % reshape to ensure horizontal vector in empty case
	EsNotMapped=reshape(EsNotMapped,1,length(EsNotMapped));
	GTsNotMapped=reshape(GTsNotMapped,1,length(GTsNotMapped));

    alldist=Inf*ones(Ngt,N);

    mindist=0;
    for o=GTsNotMapped
        GT=[gtInfo.X(t,o) gtInfo.Y(t,o)];
        rd = gtInfo.rd(t,o);
        for e=EsNotMapped
            E=[stateInfo.X(t,e) stateInfo.Y(t,e)];
            alldist(o,e)=norm(GT-E)/rd;     % divide by the reference distance
        end
    end

    tmpai=alldist;        
    tmpai(tmpai>td)=Inf; %td -> threshold for matching
    [Mtch,Cst]=Hungarian(tmpai);
    [u,v]=find(Mtch);

    for mmm=1:length(u)
        M(t,u(mmm))=v(mmm);
    end
        
    curtracked=find(M(t,:));
    
    
    alltrackers=find(stInd(t,:));
    mappedtrackers=intersect(M(t,find(M(t,:))),alltrackers);
    falsepositives=setdiff(alltrackers,mappedtrackers);
    
    alltracked(t,:)=M(t,:);
%     allfalsepos(t,1:length(falsepositives))=falsepositives;
    allfalsepos(t,falsepositives)=falsepositives;
    
    %%  mismatch errors
    if t>1
        for ct=curtracked
            lastnotempty=find(M(1:t-1,ct),1,'last');
%             if gtInd(t-1,ct) && ~isempty(lastnotempty) && M(t,ct)~=M(lastnotempty,ct)
            if ~isempty(lastnotempty) && M(t,ct)~=M(lastnotempty,ct)
                mme(t)=mme(t)+1;
            end
        end
    end
    
    c(t)=numel(curtracked);
    for ct=curtracked
        eid=M(t,ct);
        d(t,ct)=norm([gtInfo.X(t,ct) gtInfo.Y(t,ct)] - ...
                [stateInfo.X(t,eid) stateInfo.Y(t,eid)])/gtInfo.rd(t,ct);
    end
    
    
    fp(t)=numel(find(stInd(t,:)))-c(t);
    m(t)=g(t)-c(t);
end    

missed=sum(m);
falsepositives=sum(fp);
idswitches=sum(mme);
numGT = sum(g);
numC = sum(c);
numD = sum(sum(d));

% MOTP=(1-sum(sum(d))/sum(c)/td) * 100; % avg distance to [0,100] 
% 
% MOTAL=(1-((sum(m)+sum(fp)+log10(sum(mme)+1))/sum(g)))*100;
% MOTA=(1-((sum(m)+sum(fp)+(sum(mme)))/sum(g)))*100;
% recall=sum(c)/sum(g)*100;
% precision=sum(c)/(sum(fp)+sum(c))*100;
% FAR=sum(fp)/Fgt;
 

%% MT PT ML
MTstatsa=zeros(1,Ngt);
for i=1:Ngt
    gtframes=find(gtInd(:,i));
    gtlength=length(gtframes);
    gttotallength=numel(find(gtInd(:,i)));
    trlengtha=numel(find(alltracked(gtframes,i)>0));
    if gtlength/gttotallength >= 0.8 && trlengtha/gttotallength < 0.2
        MTstatsa(i)=3;
    elseif t>=find(gtInd(:,i),1,'last') && trlengtha/gttotallength <= 0.8
        MTstatsa(i)=2;
    elseif trlengtha/gttotallength >= 0.8
        MTstatsa(i)=1;
    end
end
% MTstatsa
MT=numel(find(MTstatsa==1));PT=numel(find(MTstatsa==2));ML=numel(find(MTstatsa==3));

%% fragments
fr=zeros(1,Ngt);
for i=1:Ngt
%      b=alltracked(find(alltracked(:,i),1,'first'):find(alltracked(:,i),1,'last'),i);
%      b(~~b)=1;
%      fr(i)=numel(find(diff(b)==-1));
    
    % determine GT snippets
    frags=~~gtInfo.X(:,i);
    starts=find(frags(1:end-1)==frags(2:end)-1)+1;
    ends=find(frags(1:end-1)==frags(2:end)+1);
    if frags(1), starts=[1; starts]; end
    if frags(end), ends=[ends; numel(frags)]; end
    
    % only count fragments within snippets
    for s=1:numel(starts)
      fragframes=starts(s):ends(s);
      b=alltracked(fragframes,i);
      b(~~b)=1;
      fr(i)=fr(i)+numel(find(diff(b)==-1));
    end
    
end
FRA=sum(fr);

assert(Ngt==MT+PT+ML,'Hmm... Not all tracks classified correctly.');
metrics=[Ngt, MT, PT, ML, falsepositives, missed, idswitches, FRA, numGT, numC, numD, Fgt];

additionalInfo.alltracked=alltracked;
additionalInfo.allfalsepos=allfalsepos;
end 

function ret=matched(gtInfo,stateInfo,t,map,mID,td)
    Xgt=gtInfo.X(t,map); Ygt=gtInfo.Y(t,map);
    rd=gtInfo.rd(t,map);
    X=stateInfo.X(t,mID); Y=stateInfo.Y(t,mID);
    ret=norm([Xgt Ygt]-[X Y])/rd<=td;
end