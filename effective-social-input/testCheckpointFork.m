function testCheckpointFork()
% A resumed branch must preserve the complete neural-field coordinate state.
base=defaultEffectiveInputConfig();
base.NAgents=8;
base.NNeurons=40;
base.NSteps=60;
base.burnIn=20;
base.recordStride=1;
base.seed=136;
base.scheduledOcclusion=true;
base.selectionPolicy='fixed-id';
base.neighborCount=4;
base.occlusionHiddenCount=2;
base.occlusionStartStep=2;
base.occlusionDurationSteps=0;
base.predictionMode='none';
checkpoint=simulateEffectiveSocialInput(base);

branch=base;
branch.NSteps=50;
branch.burnIn=0;
branch.initialXY=checkpoint.finalPosition;
branch.initialHeadings=checkpoint.finalHeadings;
branch.initialNeuralState=checkpoint.finalNeuralState;
branch.initialPreferredDirections=checkpoint.finalPreferredDirections;
branch.occlusionStartStep=11;
branch.occlusionDurationSteps=0;

none0=simulateEffectiveSocialInput(branch);
branch.predictionMode='hold';
hold0=simulateEffectiveSocialInput(branch);
assert(isequaln(none0.timeSeries,hold0.timeSeries));
assert(isequaln(none0.trajectory,hold0.trajectory));
assert(isequaln(none0.finalNeuralState,hold0.finalNeuralState));
assert(isequaln(none0.finalPreferredDirections,hold0.finalPreferredDirections));

branch.occlusionDurationSteps=20;
holdT=simulateEffectiveSocialInput(branch);
branch.predictionMode='none';
deleteT=simulateEffectiveSocialInput(branch);
pre=1:10;
fields=fieldnames(holdT.timeSeries);
for i=1:numel(fields)
    a=holdT.timeSeries.(fields{i})(pre);
    b=deleteT.timeSeries.(fields{i})(pre);
    assert(isequaln(a,b));
end
assert(isequaln(holdT.trajectory.x(:,1:11),deleteT.trajectory.x(:,1:11)));
assert(isequaln(holdT.trajectory.y(:,1:11),deleteT.trajectory.y(:,1:11)));
assert(isequaln(holdT.trajectory.headings(:,1:11),deleteT.trajectory.headings(:,1:11)));

repeat=simulateEffectiveSocialInput(branch);
assert(isequaln(deleteT.timeSeries,repeat.timeSeries));
assert(isequaln(deleteT.trajectory,repeat.trajectory));
disp('PASS: complete checkpoint, T=0 equality, common pre-event path, deterministic replay');
end
