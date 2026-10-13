function [mask,diagnostics]=makeProtocolMasks(position,cfg,bearingOverride)
%MAKEPROTOCOLMASKS Enumerate sorted actual-ID tuples before selecting a tie.
validateattributes(cfg.h,{'numeric'},{'scalar','integer','>=',0,'<=',cfg.k});
n=cfg.NAgents; h=cfg.h; targets=fixedTargetIds(cfg);
mask=zeros(n,h); diagnostics=cell(n,1);
stream=RandStream('mt19937ar','Seed',100000+cfg.seed);
for i=1:n
    ids=sort(targets(i,:));
    if h==0, tuples=zeros(1,0); else, tuples=nchoosek(ids,h); end
    if h>0, tuples=sortrows(tuples); end
    dxy=position(:,ids)-position(:,i);
    dxy=dxy-cfg.arenaSize*round(dxy/cfg.arenaSize);
    b=mod(atan2(dxy(2,:),dxy(1,:)),2*pi);
    if nargin==3, b=bearingOverride(i,ids); end
    R1=zeros(size(tuples,1),1); R2=R1; gaps=nan(size(R1));
    for q=1:size(tuples,1)
        [~,ix]=ismember(tuples(q,:),ids); s=angularSetMetrics(b(ix));
        R1(q)=s.R1; R2(q)=s.R2; gaps(q)=s.maxGap;
    end
    switch cfg.maskRule
        case 'random', pick=randi(stream,size(tuples,1));
        case 'max_hidden_R1', pick=find(abs(R1-max(R1))<=1e-12,1);
        case 'min_hidden_R1', pick=find(abs(R1-min(R1))<=1e-12,1);
        otherwise, error('SQ1:MaskRule','Unknown mask rule.');
    end
    mask(i,:)=tuples(pick,:);
    [~,hiddenIx]=ismember(mask(i,:),ids); remaining=~ismember(ids,mask(i,:));
    d.candidateIds=tuples; d.R1=R1; d.R2=R2; d.maxGap=gaps;
    d.chosenIndex=pick; d.chosenIds=mask(i,:); d.targetIds=targets(i,:);
    d.hidden=angularSetMetrics(b(hiddenIx));
    d.remaining=angularSetMetrics(b(remaining)); d.all=angularSetMetrics(b);
    d.generator='mt19937ar'; d.seed=100000+cfg.seed;
    diagnostics{i}=d;
end
end
