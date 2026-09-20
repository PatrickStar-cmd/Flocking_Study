function runHighRateBearingValidation()
% Four reference runs, followed by offline masked-observation forecasting.
root=fileparts(mfilename('fullpath'));
out=fullfile(root,'results','bearing-high-rate-v1');
if ~exist(out,'dir'), mkdir(out); end
seeds=138:141; qs=[1e-8 1e-6 1e-4]; R=1e-4;
starts=[6001 6801 7601]; horizons=[10 34 100 200 400];
base=defaultEffectiveInputConfig(); base.recordStride=1;
meta=struct('seeds',seeds,'trainSeeds',138:139,'evaluationSeeds',140:141, ...
    'qCandidates',qs,'R',R,'starts',starts,'horizonSteps',horizons, ...
    'dt',base.dt,'historySteps',200,'matlabVersion',version,'base',base);
metaPath=fullfile(out,'metadata.mat');
if isfile(metaPath)
    prior=load(metaPath,'meta'); assert(isequaln(prior.meta,meta));
else
    save(metaPath,'meta');
end
for seed=seeds
    path=fullfile(out,sprintf('reference-%d.mat',seed));
    if isfile(path)
        old=load(path,'complete'); assert(old.complete); continue;
    end
    cfg=base; cfg.seed=seed;
    fprintf('Generating reference seed %d\n',seed);
    result=simulateEffectiveSocialInput(cfg); complete=true;
    save([path '.partial'],'result','complete','-v7'); movefile([path '.partial'],path);
    fprintf('Saved reference seed %d\n',seed);
end
rows=cell(numel(seeds)*20*numel(starts)*numel(qs)*numel(horizons),12); idx=0;
for seed=seeds
    x=load(fullfile(out,sprintf('reference-%d.mat',seed)),'result');
    p=x.result.trajectory; dt=base.dt;
    assert(numel(p.time)==8001 && max(abs(diff(p.time)-dt))<1e-10);
    stream=RandStream('mt19937ar','Seed',seed+9000);
    for i=1:20
        j=mod(i,20)+1; L=base.arenaSize;
        dx=p.x(j,:)-p.x(i,:); dy=p.y(j,:)-p.y(i,:);
        dx=dx-L*round(dx/L); dy=dy-L*round(dy/L);
        truth=atan2(dy,dx);
        observation=truth+sqrt(R)*randn(stream,size(truth));
        jump=[false,abs(diff(dx))>L/2 | abs(diff(dy))>L/2];
        for start=starts
            for q=qs
                track=[]; nis=[];
                for ti=start-199:start
                    track=bearingTrackStep(track,observation(ti),dt,q,R);
                    if ti>start-100, nis(end+1)=track.nis; end %#ok<AGROW>
                end
                for h=1:max(horizons)
                    track=bearingTrackStep(track,[],dt,q,R);
                    if ~ismember(h,horizons), continue; end
                    ek=atan2(sin(track.x(1)-truth(start+h)),cos(track.x(1)-truth(start+h)));
                    eh=atan2(sin(observation(start)-truth(start+h)),cos(observation(start)-truth(start+h)));
                    variance=track.P(1,1); idx=idx+1;
                    rows(idx,:)={seed,i,start,q,h*dt,ek^2,eh^2,variance, ...
                        abs(ek)<=1.96*sqrt(variance),1.96*sqrt(variance)>=pi, ...
                        any(jump(start+1:start+h)),mean(nis)};
                end
            end
        end
    end
    fprintf('Scored seed %d\n',seed);
end
t=cell2table(rows,'VariableNames',{'seed','focal','startIndex','q','T','kalmanSquaredError', ...
    'holdSquaredError','variance','covered','intervalCoversCircle','periodicJump','historyNIS'});
writetable(t,fullfile(out,'errors.csv'));
loss=zeros(size(qs));
for qi=1:numel(qs), loss(qi)=mean(t.kalmanSquaredError(t.seed<=139 & t.q==qs(qi))); end
[~,qi]=min(loss); selectedQ=qs(qi); reports=cell(5,9);
for hi=1:numel(horizons)
    v=t(t.seed>=140 & t.q==selectedQ & t.T==horizons(hi)*base.dt,:);
    improvements=zeros(2,1);
    for si=1:2
        mask=v.seed==139+si;
        improvements(si)=mean(v.kalmanSquaredError(mask)-v.holdSquaredError(mask));
    end
    reports(hi,:)={horizons(hi)*base.dt,sqrt(mean(v.kalmanSquaredError)),sqrt(mean(v.holdSquaredError)), ...
        mean(v.covered),mean(v.intervalCoversCircle),sum(v.periodicJump),sum(improvements<0),improvements(1),improvements(2)};
end
summary=cell2table(reports,'VariableNames',{'T','kalmanRMSE','holdRMSE','coverage95', ...
    'wholeCircleFraction','jumpCasesOf120','improvedSeedsOf2','seed140MSEDifference','seed141MSEDifference'});
writetable(summary,fullfile(out,'summary.csv'));
meta.selectedQ=selectedQ; meta.trainingMSE=loss;
fid=fopen(fullfile(out,'selection.json'),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
disp(selectedQ); disp(loss); disp(summary);
end
