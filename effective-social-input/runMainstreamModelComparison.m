function runMainstreamModelComparison()
%RUNMAINSTREAMMODELCOMPARISON Representative mainstream flocking baselines.
% All methods start from the same neural-model checkpoint and experience the
% same directed observation blackout. Baselines use richer state information
% and are therefore reference comparators, not bearing-only replacements.
root=fileparts(mfilename('fullpath'));
kfDir=fullfile(root,'results','paper-kf-simulation');
out=fullfile(root,'results','mainstream-model-comparison');
if ~exist(out,'dir'), mkdir(out); end
seeds=101:105; models={'allocentric-kf-unc','vicsek','boids','olfati-saber-type'};
Tsteps=100; warm=200; recover=300; total=warm+Tsteps+recover; event=(warm+1):(warm+Tsteps); post=(warm+Tsteps+1):total;
rows={}; idx=0;
for seed=seeds
 cp=load(fullfile(kfDir,sprintf('checkpoint-%d.mat',seed)),'checkpoint'); cp=cp.checkpoint;
 speed=mean(cp.timeSeries.meanSpeed(end-199:end));
 % Proposed method, rerun with and without blackout from the same checkpoint.
 cfg=cp.config; cfg.NSteps=total; cfg.burnIn=0; cfg.recordStride=1;
 cfg.initialXY=cp.finalPosition; cfg.initialHeadings=cp.finalHeadings;
 cfg.initialNeuralState=cp.finalNeuralState; cfg.initialPreferredDirections=cp.finalPreferredDirections;
 cfg.occlusionStartStep=warm+1; cfg.predictionMode='kalman';
 cfg.occlusionDurationSteps=0; c=simulateEffectiveSocialInput(cfg);
 cfg.occlusionDurationSteps=Tsteps; r=simulateEffectiveSocialInput(cfg);
 idx=idx+1; rows(idx,:)=makeNeuralRow(seed,models{1},c,r,event,post);
 for mi=2:numel(models)
   c=simulateRuleModel(cp.finalPosition,cp.finalHeadings,speed,cfg,models{mi},warm,0,total,seed);
   r=simulateRuleModel(cp.finalPosition,cp.finalHeadings,speed,cfg,models{mi},warm,Tsteps,total,seed);
   idx=idx+1; rows(idx,:)=makeRuleRow(seed,models{mi},c,r,event,post);
 end
 fprintf('mainstream comparison seed %d complete\n',seed);
end
R=cell2table(rows,'VariableNames',{'seed','model','eventGO','eventGOLoss','postPairDistance','postPairDistanceIncrease','minimumDistance','collisionFraction','controlEffort','recoveryTime','informationLevel'});
writetable(R,fullfile(out,'summary.csv'));
fid=fopen(fullfile(out,'README.md'),'w');
fprintf(fid,['# Mainstream model comparison\n\nPaired exploratory comparison, seeds 101--105, N=20, T=30. ' ...
 'All models start from the same allocentric-neural-field checkpoint. Vicsek uses neighbor headings; Boids and the Olfati--Saber-type potential-consensus baseline use relative positions and headings. ' ...
 'They therefore have richer observations than the proposed bearing-only method. Parameters are representative and not independently optimized; conclusions should emphasize sensing-efficiency and robustness patterns, not a universal algorithm ranking.\n']); fclose(fid);
end

function row=makeNeuralRow(seed,name,c,r,event,post)
go=r.timeSeries.globalOrder; cg=c.timeSeries.globalOrder;
pd=r.timeSeries.meanPairDistanceNorm; cpd=c.timeSeries.meanPairDistanceNorm;
 [minD,coll]=trajectorySafetyStats(r.trajectory,r.config.arenaSize,10);
eff=mean(abs(diff(r.trajectory.headings,1,2)),'all');
rec=recoveryTime(go,cg,post,r.config.dt);
 row={seed,name,mean(go(event)),mean(cg(event)-go(event)),mean(pd(post)),mean(pd(post)-cpd(post)),minD,coll,eff,rec,'bearing only'};
end

function row=makeRuleRow(seed,name,c,r,event,post)
rec=recoveryTime(r.go,c.go,post,r.dt);
row={seed,name,mean(r.go(event)),mean(c.go(event)-r.go(event)),mean(r.pair(post)),mean(r.pair(post)-c.pair(post)),min(r.minimum)/r.L,mean(r.minimum<10),mean(r.effort),rec,r.information};
end

function rec=recoveryTime(x,control,post,dt)
tol=0.05; ok=abs(x(post)-control(post))<=tol; j=find(conv(double(ok),ones(10,1),'same')>=10,1);
if isempty(j), rec=NaN; else, rec=(j-1)*dt; end
end

function out=simulateRuleModel(pos,head,speed,cfg,model,warm,duration,total,seed)
rng(seed+sum(double(model)),'twister'); n=size(pos,2); L=cfg.arenaSize; dt=cfg.dt;
go=zeros(total,1); pair=zeros(total,1); minimum=zeros(total,1); effort=zeros(total,1);
for step=1:total
 old=head; desired=head;
 hidden=step>warm && step<=warm+duration;
 for i=1:n
   ids=mod(i-1+(1:cfg.neighborCount),n)+1;
   if hidden, ids=ids((cfg.occlusionHiddenCount+1):end); end
   dx=pos(:,ids)-pos(:,i); dx=dx-L*round(dx/L); dist=hypot(dx(1,:),dx(2,:));
   switch model
    case 'vicsek'
     z=exp(1i*head(ids)); desired(i)=angle(exp(1i*head(i))+sum(z));
    case 'boids'
     align=[mean(cos(head(ids)));mean(sin(head(ids)))];
     cohesion=mean(dx,2)/L; near=dist<30 & dist>0;
     separation=[0;0]; if any(near), separation=-sum(dx(:,near)./(dist(near).^2+eps),2); end
     v=1.0*align+2.0*cohesion+15*separation; desired(i)=atan2(v(2),v(1));
    otherwise
     align=[mean(cos(head(ids)));mean(sin(head(ids)))];
     d0=80; potential=sum(dx.*((dist-d0)./(dist+eps)),2)/max(1,numel(ids))/d0;
     repel=sum(-dx./(dist.^2+eps),2);
     v=align+0.7*potential+20*repel; desired(i)=atan2(v(2),v(1));
   end
 end
 gain=0.25; delta=atan2(sin(desired-head),cos(desired-head)); head=mod(head+gain*delta,2*pi);
 pos=mod(pos+dt*speed*[cos(head)';sin(head)'],L);
 effort(step)=mean(abs(atan2(sin(head-old),cos(head-old))))/dt;
 go(step)=abs(mean(exp(1i*head)));
 [pair(step),minimum(step)]=distanceStats(pos,L);
end
out=struct('go',go,'pair',pair/L,'minimum',minimum,'effort',effort,'dt',dt,'L',L);
if strcmp(model,'vicsek'), out.information='neighbor heading'; else, out.information='relative position + heading'; end
end

function [mp,mn]=distanceStats(pos,L)
n=size(pos,2); vals=[]; nearest=inf(n,1);
for i=1:n-1
 d=pos(:,i+1:n)-pos(:,i); d=d-L*round(d/L); z=hypot(d(1,:),d(2,:)); vals=[vals z];
 nearest(i)=min(nearest(i),min(z));
 for j=i+1:n, nearest(j)=min(nearest(j),z(j-i)); end
end

mp=mean(vals); mn=min(nearest);
end

function [minNorm,collisionFraction]=trajectorySafetyStats(tr,L,threshold)
minD=inf; collisions=0; count=0;
for k=1:size(tr.x,2)
 p=[tr.x(:,k)';tr.y(:,k)']; [~,d]=distanceStats(p,L); minD=min(minD,d); collisions=collisions+(d<threshold); count=count+1;
end
minNorm=minD/L; collisionFraction=collisions/max(1,count);
end
