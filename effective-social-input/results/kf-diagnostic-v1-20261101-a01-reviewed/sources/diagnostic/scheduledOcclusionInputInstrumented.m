function [input,state,adj,m,record] = scheduledOcclusionInputInstrumented(position,alpha,state,cfg,step)
%SCHEDULEDOCCLUSIONINPUTINSTRUMENTED Mechanically identical to stock
% scheduledOcclusionInput, with added PASSIVE diagnostic recording.
% Hidden truth is NEVER passed to the controller.

if isempty(state)
    stream = RandStream('mt19937ar', 'Seed', cfg.seed);
    sensorNoise = cfg.bearingMeasurementStd * randn(stream, cfg.NAgents, cfg.NAgents, cfg.NSteps);
else
    sensorNoise = state.sensorNoise;
end

n = cfg.NAgents;
selected = false(n);
fresh = false(n);
bearings = nan(n);
weights = nan(n);
trueBearings = nan(n);
hidden = step >= cfg.occlusionStartStep && step < cfg.occlusionStartStep + cfg.occlusionDurationSteps;

for i = 1:n
    ids = mod(i - 1 + (1:cfg.neighborCount), n) + 1;
    selected(i, ids) = true;
    count = min(cfg.occlusionHiddenCount, numel(ids));
    visible = ids;
    if hidden
        if isempty(cfg.occlusionHiddenIds)
            visible = ids(count + 1:end);
        else
            hiddenIds = cfg.occlusionHiddenIds(i, 1:count);
            visible = ids(~ismember(ids, hiddenIds));
        end
    end
    fresh(i, visible) = true;
    dx = position(:, visible) - position(:, i);
    dx = dx - cfg.arenaSize * round(dx / cfg.arenaSize);
    bearings(i, visible) = mod(atan2(dx(2, :), dx(1, :)) + sensorNoise(i, visible, step), 2*pi);
    weights(i, visible) = distanceWeight(hypot(dx(1, :), dx(2, :)), cfg);
    dxAll = position(:, ids) - position(:, i);
    dxAll = dxAll - cfg.arenaSize * round(dxAll / cfg.arenaSize);
    trueBearings(i, ids) = mod(atan2(dxAll(2, :), dxAll(1, :)), 2*pi);
end

obs = struct('selected', selected, 'fresh', fresh, 'bearing', bearings, 'weight', weights);
opt = struct('dt', cfg.dt, 'q', cfg.bearingProcessNoise, ...
    'gmQ', cfg.bearingGaussMarkovProcessNoise, 'velocityTau', cfg.bearingVelocityTau, ...
    'R', cfg.bearingObservationVariance, ...
    'mode', cfg.predictionMode, 'maxAge', cfg.predictionMaxAge, 'maxVariance', cfg.predictionMaxVariance, ...
    'sigma', cfg.receptiveFieldWidth, 'h0', cfg.totalSocialAttraction / n);

if ~isempty(state), state = rmfield(state, 'sensorNoise'); end
[input, state, d] = occlusionInputStepInstrumented(obs, alpha, state, opt);
adj = d.active;
state.sensorNoise = sensorNoise;

m.visibleCount = sum(fresh, 2);
m.activeCount = sum(adj, 2);
m.amplitudeSum = sum(adj .* fillmissing(state.weight, 'constant', 0), 2) * opt.h0;
m.inputMass = sum(input, 1)' * 2 * pi / cfg.NNeurons;
m.fieldPeak = max(input, [], 1)';
m.distanceWeight = sum(adj .* fillmissing(state.weight, 'constant', 0), 2) ./ max(1, m.activeCount);
z = sum(input .* exp(1i * alpha'), 1);
m.inputConcentration = abs(z)' ./ max(eps, sum(input, 1)');
m.angularCoverage = nan(n, 1);

% Stock record fields (100% backward compatible):
record = struct('fresh', fresh, 'active', logical(adj), 'variance', d.variance, 'age', d.age, 'weight', d.weight);
record.observedBearing = bearings;
record.injectedBearing = d.center;
record.trueBearing = trueBearings;

% Passive diagnostic fields:
record.rawVariance = d.rawVariance;
record.rawP22 = d.rawP22;
record.rawP12 = d.rawP12;
record.rawP21 = d.rawP21;
record.rawX = d.rawX;
record.rawAge = d.rawAge;
record.rawInnovation = d.rawInnovation;
record.rawNis = d.rawNis;
record.cutoffReason = d.cutoffReason;
end
