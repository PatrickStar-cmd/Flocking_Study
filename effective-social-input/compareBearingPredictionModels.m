function compareBearingPredictionModels()
% Offline model comparison using saved high-rate reference trajectories.
root=fileparts(mfilename('fullpath'));
sourceDir=fullfile(root,'results','bearing-high-rate-v1');
out=fullfile(root,'results','bearing-damped-model-v1');
if ~exist(out,'dir'), mkdir(out); end

seeds=138:141;
trainSeeds=138:139;
evaluationSeeds=140:141;
qs=[1e-6 1e-5 1e-4 1e-3];
taus=[1 3 10 30 60 Inf];
R=1e-4;
starts=[6001 6801 7601];
horizons=[10 34 100 200 400];
dt=.3;

rows=cell(numel(seeds)*20*numel(starts)*numel(qs)*numel(taus)*numel(horizons),14);
idx=0;
for seed=seeds
    loaded=load(fullfile(sourceDir,sprintf('reference-%d.mat',seed)),'result');
    p=loaded.result.trajectory;
    assert(numel(p.time)==8001 && max(abs(diff(p.time)-dt))<1e-10);
    stream=RandStream('mt19937ar','Seed',seed+9000);
    for i=1:20
        j=mod(i,20)+1;
        L=loaded.result.config.arenaSize;
        dx=p.x(j,:)-p.x(i,:);
        dy=p.y(j,:)-p.y(i,:);
        dx=dx-L*round(dx/L);
        dy=dy-L*round(dy/L);
        truth=atan2(dy,dx);
        observation=truth+sqrt(R)*randn(stream,size(truth));
        jump=[false,abs(diff(dx))>L/2 | abs(diff(dy))>L/2];
        for start=starts
            for q=qs
                for tau=taus
                    track=[];
                    nis=[];
                    for ti=start-199:start
                        track=bearingTrackStepGaussMarkov(track,observation(ti),dt,q,R,tau);
                        if ti>start-100, nis(end+1)=track.nis; end %#ok<AGROW>
                    end
                    startRate=track.x(2);
                    for h=1:max(horizons)
                        track=bearingTrackStepGaussMarkov(track,[],dt,q,R,tau);
                        if ~ismember(h,horizons), continue; end
                        ep=atan2(sin(track.x(1)-truth(start+h)),cos(track.x(1)-truth(start+h)));
                        eh=atan2(sin(observation(start)-truth(start+h)),cos(observation(start)-truth(start+h)));
                        variance=track.P(1,1);
                        idx=idx+1;
                        rows(idx,:)={seed,i,start,q,tau,h*dt,ep^2,eh^2,variance, ...
                            abs(ep)<=1.96*sqrt(variance),1.96*sqrt(variance)>=pi, ...
                            any(jump(start+1:start+h)),mean(nis),startRate};
                    end
                end
            end
        end
    end
    fprintf('Scored seed %d\n',seed);
end
t=cell2table(rows,'VariableNames',{'seed','focal','startIndex','q','tau','T', ...
    'predictionSquaredError','holdSquaredError','variance','covered', ...
    'intervalCoversCircle','periodicJump','historyNIS','estimatedStartRate'});
writetable(t,fullfile(out,'errors.csv'));

% The explicit values avoid selecting on long horizons whose intervals span the circle.
selectionMask=ismember(t.seed,trainSeeds) & ismember(t.T,[3 10.2 30]);
selection=cell(numel(qs)*numel(taus),5);
sidx=0;
for q=qs
    for tau=taus
        v=t(selectionMask & t.q==q & t.tau==tau,:);
        sidx=sidx+1;
        selection(sidx,:)={q,tau,mean(v.predictionSquaredError),mean(v.historyNIS),height(v)};
    end
end
selection=cell2table(selection,'VariableNames',{'q','tau','trainingMSE','meanHistoryNIS','nForecasts'});
selection=sortrows(selection,{'trainingMSE','tau','q'});
writetable(selection,fullfile(out,'model-selection.csv'));

bestDamped=selection(isfinite(selection.tau),:);
bestDamped=bestDamped(1,:);
bestCV=selection(isinf(selection.tau),:);
bestCV=bestCV(1,:);
reports=cell(numel(horizons),15);
for hi=1:numel(horizons)
    T=horizons(hi)*dt;
    gm=t(ismember(t.seed,evaluationSeeds) & t.q==bestDamped.q & t.tau==bestDamped.tau & t.T==T,:);
    cv=t(ismember(t.seed,evaluationSeeds) & t.q==bestCV.q & isinf(t.tau) & t.T==T,:);
    assert(height(gm)==120 && height(cv)==120);
    diffs=zeros(2,2);
    for si=1:2
        diffs(si,1)=mean(gm.predictionSquaredError(gm.seed==evaluationSeeds(si))-gm.holdSquaredError(gm.seed==evaluationSeeds(si)));
        diffs(si,2)=mean(gm.predictionSquaredError(gm.seed==evaluationSeeds(si))-cv.predictionSquaredError(cv.seed==evaluationSeeds(si)));
    end
    reports(hi,:)={T,sqrt(mean(gm.predictionSquaredError)),sqrt(mean(cv.predictionSquaredError)), ...
        sqrt(mean(gm.holdSquaredError)),mean(gm.covered),mean(gm.intervalCoversCircle), ...
        sum(gm.periodicJump),sum(diffs(:,1)<0),sum(diffs(:,2)<0), ...
        diffs(1,1),diffs(2,1),diffs(1,2),diffs(2,2),bestDamped.q,bestDamped.tau};
end
summary=cell2table(reports,'VariableNames',{'T','dampedRMSE','cvRMSE','holdRMSE', ...
    'dampedCoverage95','dampedWholeCircleFraction','jumpCasesOf120', ...
    'dampedBetterThanHoldSeedsOf2','dampedBetterThanCVSeedsOf2', ...
    'seed140DampedMinusHoldMSE','seed141DampedMinusHoldMSE', ...
    'seed140DampedMinusCVMSE','seed141DampedMinusCVMSE','selectedQ','selectedTau'});
writetable(summary,fullfile(out,'summary.csv'));

audit=struct('createdAt',char(datetime('now','TimeZone','Asia/Shanghai')), ...
    'sourceDirectory','results/bearing-high-rate-v1','sourceSeeds',seeds, ...
    'trainSeeds',trainSeeds,'evaluationSeeds',evaluationSeeds, ...
    'qCandidates',qs,'tauCandidates',taus,'R',R,'dt',dt,'starts',starts, ...
    'horizonSteps',horizons,'selectionT',[3 10.2 30], ...
    'selectedDampedQ',bestDamped.q,'selectedDampedTau',bestDamped.tau, ...
    'selectedCVQ',bestCV.q,'matlabVersion',version);
fid=fopen(fullfile(out,'audit.json'),'w');
fprintf(fid,'%s',jsonencode(audit,'PrettyPrint',true));
fclose(fid);
disp(selection(1:min(10,height(selection)),:));
disp(summary);
end
