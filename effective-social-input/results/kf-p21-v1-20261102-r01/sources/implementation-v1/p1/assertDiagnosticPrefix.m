function assertDiagnosticPrefix(reference,other,lastPreStep)
% Pre-state record t is state after t-1; include t0 to check last fresh update.
idx=1:(lastPreStep+1);
assert(isequaln(reference.positionsRecord(:,:,idx),other.positionsRecord(:,:,idx)));
assert(isequaln(reference.headingsRecord(:,idx),other.headingsRecord(:,idx)));
assert(isequaln(reference.uRecord(:,:,idx),other.uRecord(:,:,idx)));
assert(isequaln(reference.finalPreferredDirections,other.finalPreferredDirections));
assert(isequaln(reference.sensorNoise,other.sensorNoise));
assert(isequaln(reference.preOcclusionCheckpoint,other.preOcclusionCheckpoint));
names=fieldnames(reference.timeSeries);
for i=1:numel(names)
    key=names{i};
    assert(isequaln(reference.timeSeries.(key)(1:lastPreStep), ...
        other.timeSeries.(key)(1:lastPreStep)),'Prefix metric mismatch %s',key);
end
assert(reference.preOcclusionCheckpoint.stateAfterStep==lastPreStep);
end
