function runDistanceWeightPilot()
root=fileparts(mfilename('fullpath')); out=fullfile(root,'results','distance-weight-pilot-v1');
if ~exist(out,'dir'), mkdir(out); end
base=defaultEffectiveInputConfig(); base.recordStride=100;
conditions={'equal','original';'exponential','original';'equal','fixed-total';'exponential','fixed-total'};
seeds=131:135; ks=[7 15]; matlabVersion=version;
save(fullfile(out,'metadata.mat'),'base','conditions','seeds','ks','matlabVersion');
for seed=seeds
 for k=ks
  for c=1:4
   path=fullfile(out,sprintf('seed-%d-k-%d-condition-%d.mat',seed,k,c));
   if isfile(path), x=load(path,'complete'); assert(x.complete); continue; end
   cfg=base; cfg.seed=seed; cfg.selectionPolicy='random'; cfg.neighborCount=k;
   cfg.distanceWeightMode=conditions{c,1}; cfg.strengthMode=conditions{c,2};
   fprintf('Start seed=%d k=%d condition=%d\n',seed,k,c);
   result=simulateEffectiveSocialInput(cfg); complete=true;
   save([path '.partial'],'result','complete','-v7'); movefile([path '.partial'],path);
   fprintf('Saved seed=%d k=%d condition=%d\n',seed,k,c);
  end
 end
end
disp('COMPLETE: 40 distance-weight pilot simulations');
end
