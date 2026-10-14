function validate_reference(out)
% Reference only: calls the original repository APIs without modifying them.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
request=jsondecode(fileread(fullfile(out,'request.json')));
c=defaultEffectiveInputConfig();
c.NSteps=12; c.burnIn=0; c.recordStride=1; c.NAgents=20;
c.NNeurons=100; c.neighborCount=15; c.selectionPolicy='fixed-id';
c.scheduledOcclusion=true; c.occlusionDurationSteps=0;
c.bearingMeasurementStd=0; c.initialXY=request.xy;
c.initialHeadings=request.heading; c.initialPreferredDirections=request.alpha;
c.initialNeuralState=request.u; c.predictionMode='none';
r=simulateEffectiveSocialInput(c);
reference=struct('u',r.finalNeuralState,'x',r.trajectory.x, ...
 'y',r.trajectory.y,'heading',r.trajectory.headings);
samples=[6.25,6.27,.01,.025,.04,.055,nan(1,8),.19];
reference.samples=samples;
alpha=repmat((0:39)*2*pi/40,2,1);
opt=struct('dt',.3,'q',1e-4,'gmQ',1e-5,'velocityTau',10, ...
 'R',1e-4,'maxAge',120,'maxVariance',1,'sigma',.4,'h0',.012);
for mode={'hold','mean','kalman'}
 opt.mode=mode{1}; state=[]; values=cell(1,numel(samples));
 for j=1:numel(samples)
  obs=struct('selected',[false true;false false], ...
   'fresh',false(2),'bearing',nan(2),'weight',nan(2));
  if isfinite(samples(j))
   obs.fresh(1,2)=true; obs.bearing(1,2)=samples(j); obs.weight(1,2)=1;
  end
  [field,state]=occlusionInputStep(obs,alpha,state,opt);
  values{j}=struct('field',field,'x',state.tracks{1,2}.x,'P',state.tracks{1,2}.P);
 end
 reference.filters.(mode{1})=values;
end
f=fopen(fullfile(out,'matlab-reference.json'),'w');
fprintf(f,'%s',jsonencode(reference)); fclose(f);
testOcclusionPrediction;
fprintf('PASS: MATLAB reference exported; existing occlusion regression passed.\n');
end
