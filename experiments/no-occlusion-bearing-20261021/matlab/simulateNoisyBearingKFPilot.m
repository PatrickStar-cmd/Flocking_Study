function result = simulateNoisyBearingKFPilot(cfg)
%SIMULATENOISYBEARINGKFPILOT Isolated copy of the archived dynamics for a noisy-observation pilot.
% The only input-path change is noisyBearingKFPilotInputStep. Source and
% original defaults are snapshotted by the runner; protected sources stay untouched.
%
% This implementation follows the neural-field update in Salahshour and
% Couzin (2025), while exposing the number and angular selection of social
% bearing inputs. It intentionally does not add range sensing or prediction.

assert(isfield(cfg,'noisyBearingMode') && any(strcmp(cfg.noisyBearingMode,{'raw','kf-posterior'})));
assert(~cfg.scheduledOcclusion && ~cfg.turnExperiment && ~cfg.openArena);
assert(strcmp(cfg.selectionPolicy,'random') && strcmp(cfg.distanceWeightMode,'equal'));
cfg = validateConfig(cfg);
rng(cfg.seed, 'twister');

n = cfg.NAgents;
nS = cfg.NNeurons;
T = cfg.NSteps;
L = cfg.arenaSize;

if isempty(cfg.initialXY)
    positions = rand(2, n) * L;
else
    positions = cfg.initialXY;
end
if isempty(cfg.initialHeadings)
    headings = zeros(n, 1);
else
    headings = cfg.initialHeadings(:);
end
if isempty(cfg.initialNeuralState)
    u = zeros(nS, n);
else
    u = cfg.initialNeuralState;
end

% Match the authors' initialization order exactly: heading and neural
% noise are drawn agent by agent, rather than in two vectorized blocks.
if isempty(cfg.initialHeadings) || isempty(cfg.initialNeuralState)
    for agent = 1:n
        if isempty(cfg.initialHeadings)
            headings(agent) = 2*pi*rand;
        end
        if isempty(cfg.initialNeuralState)
            u(:, agent) = 0.1*randn(nS, 1);
        end
    end
end

W = buildRingConnectivity(nS, cfg.connectivityExponent);
initialPosition = positions;
initialHeadings = headings;
initialNeuralState = u;
cueTrace = nan(T,1);
alpha0 = linspace(0, 2*pi, nS + 1);
alpha0 = alpha0(1:end-1);
if isempty(cfg.initialPreferredDirections)
    alpha = mod(repmat(alpha0, n, 1) + headings, 2*pi);
else
    alpha = cfg.initialPreferredDirections;
end

recordSteps = unique([0:cfg.recordStride:T, T]);
nRecords = numel(recordSteps);
xRecord = zeros(n, nRecords);
yRecord = zeros(n, nRecords);
headingRecord = zeros(n, nRecords);
xRecord(:, 1) = positions(1, :).';
yRecord(:, 1) = positions(2, :).';
headingRecord(:, 1) = headings;
nextRecord = 2;

ts.globalOrder = nan(T, 1);
ts.localOrder = nan(T, 1);
ts.meanPairDistance = nan(T, 1);
ts.clusterRadius = nan(T, 1);
ts.meanPairDistanceNorm = nan(T, 1);
ts.meanNearestDistanceNorm = nan(T, 1);
ts.meanInputAmplitude = nan(T, 1);
ts.meanVisibleNeighbors = nan(T, 1);
ts.meanActiveNeighbors = nan(T, 1);
ts.meanAngularCoverage = nan(T, 1);
ts.meanInputConcentration = nan(T, 1);
ts.meanFieldPeak = nan(T, 1);
ts.meanBumpConcentration = nan(T, 1);
ts.meanSpeed = nan(T, 1);
ts.meanInputMass = nan(T,1);
ts.meanDistanceWeight = nan(T,1);
ts.perReceiverInputMass = nan(T, n);
ts.perReceiverHiddenWeight = nan(T, n);
ts.inputSubstitutionRms=nan(T,1);
inputSamples={};
if cfg.scheduledOcclusion
    % Prediction-quality diagnostics: computed from the injected field and the
    % evaluator's hidden truth. They never feed back into the controller.
    ts.meanPredictionError = nan(T,1);
    ts.meanPredictionSquaredError = nan(T,1);
    ts.meanPredictionVariance = nan(T,1);
end
if cfg.storeInputDiagnostics
    ts.meanInputMassExact = nan(T,1);
    ts.meanInputHarmonic1 = nan(T,1);
    ts.meanInputDirectionalInfo = nan(T,1);
    if cfg.scheduledOcclusion
        ts.meanInputMassRelError = nan(T,1);
    end
end

ts.observedBearingRMSE = nan(T,1);
ts.usedBearingRMSE = nan(T,1);
ts.inputTruthRMS = nan(T,1);
ts.meanKFVariance = nan(T,1);
pilotState = [];
adjacencySum = zeros(n, n);
evaluationCount = 0;
cosAlpha = cos(alpha).';
sinAlpha = sin(alpha).';

for t = 1:T
    sensorHeadings = headings;
    if cfg.scheduledOcclusion
        if t==1
            % Carry the sensor/filter state across a fork when one is given.
            % This is what makes the sensor-noise realisation shared: the
            % branch keeps the SAME noise array the prefix generated, so step
            % k of the branch reads the noise of absolute step prefix+k.
            if isfield(cfg,'observationState') && ~isempty(cfg.observationState)
                observationState = cfg.observationState;
            else
                observationState = [];
            end
            if cfg.storeObservationRecords, observationRecords=cell(T,1); end
        end
        if cfg.storeObservationRecords
            [externalInput,observationState,adjacency,inputMetrics,observationRecords{t}]= ...
                scheduledOcclusionInput(positions,alpha,observationState,cfg,t);
        else
            [externalInput,observationState,adjacency,inputMetrics]= ...
                scheduledOcclusionInput(positions,alpha,observationState,cfg,t);
        end
    else
        [externalInput, adjacency, inputMetrics, pilotState] = noisyBearingKFPilotInputStep( ...
            positions, alpha, pilotState, cfg, t);
        ts.observedBearingRMSE(t)=inputMetrics.observedBearingRMSE;
        ts.usedBearingRMSE(t)=inputMetrics.usedBearingRMSE;
        ts.inputTruthRMS(t)=inputMetrics.inputTruthRMS;
        ts.meanKFVariance(t)=inputMetrics.meanKFVariance;
    end
    if cfg.scheduledOcclusion
        ts.inputSubstitutionRms(t)=inputMetrics.inputSubstitutionRms;
        if cfg.storeInputDiagnostics && mod(t-1,10)==0
            inputSamples{end+1}=struct('step',t,'actual',inputMetrics.actualInput, ...
                'truth',inputMetrics.truthInput); %#ok<AGROW>
        end
    end

    if cfg.turnExperiment && cfg.turnCueAmplitude > 0 && ...
            cfg.absoluteStepOffset+t >= cfg.turnCueStartStep && ...
            cfg.absoluteStepOffset+t < cfg.turnCueStartStep + cfg.turnCueDurationSteps
        % absoluteStepOffset lets a forked run continue an intervention whose
        % clock started in the shared prefix, so the cue phase is continuous
        % across the fork instead of restarting.
        % task2-relative-confidence-v1: support in-branch piecewise cue rate schedule
        if isfield(cfg, 'turnCueSchedule') && ~isempty(cfg.turnCueSchedule) && ...
                t >= cfg.turnCueSchedule.splitStep
            tSplit = cfg.turnCueSchedule.splitStep;
            elapsedSplit = (cfg.absoluteStepOffset + (tSplit - 1) - cfg.turnCueStartStep + 1) * cfg.dt;
            elapsedAfter = (t - (tSplit - 1)) * cfg.dt;
            for leader = cfg.turnLeaderIds(:).'
                baseSplit = cfg.turnCueBaseAngles(leader) + cfg.turnCueRate * elapsedSplit;
                target = mod(baseSplit + cfg.turnCueSchedule.rate2 * elapsedAfter, 2*pi);
                cueTrace(t) = target;
                delta = atan2(sin(alpha(leader,:)-target), cos(alpha(leader,:)-target));
                externalInput(:,leader) = externalInput(:,leader) + ...
                    cfg.turnCueAmplitude*exp(-0.5*(delta/cfg.receptiveFieldWidth).^2).';
            end
        else
            elapsed = (cfg.absoluteStepOffset + t - cfg.turnCueStartStep + 1) * cfg.dt;
            for leader = cfg.turnLeaderIds(:).'
                target = mod(cfg.turnCueBaseAngles(leader) + cfg.turnCueRate*elapsed, 2*pi);
                cueTrace(t) = target;
                delta = atan2(sin(alpha(leader,:)-target), cos(alpha(leader,:)-target));
                externalInput(:,leader) = externalInput(:,leader) + ...
                    cfg.turnCueAmplitude*exp(-0.5*(delta/cfg.receptiveFieldWidth).^2).';
            end
        end
    end

    % NOTE on stochasticity: the neural field is NOT perturbed per step. The
    % authors' model draws neural noise only at initialization (agent by
    % agent), so a run is deterministic given its seed, initial state, and
    % observation stream. Task 3 therefore needs no per-step neural noise
    % stream: the shared-prefix guarantee comes entirely from the
    % pre-generated sensor noise indexed by absolute step.
    neuralOutput = tanh(cfg.beta*u);
    u = u + cfg.dt*(-u + (W*neuralOutput)/nS ...
        - cfg.globalInhibition + externalInput);

    positiveOutput = tanh(cfg.beta*u);
    positiveOutput(positiveOutput < 0) = 0;
    cx = sum(positiveOutput .* cosAlpha, 1);
    cy = sum(positiveOutput .* sinAlpha, 1);

    oldHeadings = headings;
    moving = abs(cx) + abs(cy) >= 1e-9;
    headings(moving) = mod(atan2(cy(moving), cx(moving)), 2*pi).';

    if ~cfg.allocentric
        alpha = mod(alpha - oldHeadings + headings, 2*pi);
        cosAlpha = cos(alpha).';
        sinAlpha = sin(alpha).';
        cx = sum(positiveOutput .* cosAlpha, 1);
        cy = sum(positiveOutput .* sinAlpha, 1);
    end

    velocityX = cfg.speedScale*cx;
    velocityY = cfg.speedScale*cy;
    if isfield(cfg, 'openArena') && cfg.openArena
        positions(1, :) = positions(1, :) + cfg.dt*velocityX;
        positions(2, :) = positions(2, :) + cfg.dt*velocityY;
    else
        positions(1, :) = mod(positions(1, :) + cfg.dt*velocityX, L);
        positions(2, :) = mod(positions(2, :) + cfg.dt*velocityY, L);
    end

    speed = hypot(velocityX, velocityY);
    unitX = zeros(1, n);
    unitY = zeros(1, n);
    nonzero = speed > 1e-12;
    unitX(nonzero) = velocityX(nonzero)./speed(nonzero);
    unitY(nonzero) = velocityY(nonzero)./speed(nonzero);
    ts.globalOrder(t) = hypot(sum(unitX), sum(unitY))/n;
    ts.groupHeading(t) = atan2(sum(unitY), sum(unitX));
    isOpen = isfield(cfg, 'openArena') && cfg.openArena;
    ts.localOrder(t) = topologicalLocalOrder(positions, unitX, unitY, L, 5, isOpen);

    [meanPair, meanNearest] = periodicDistanceMetrics(positions, L, isOpen);
    ts.meanPairDistance(t) = meanPair;
    ts.meanPairDistanceNorm(t) = meanPair/L;
    ts.meanNearestDistanceNorm(t) = meanNearest/L;
    ts.clusterRadius(t) = sqrt(mean(sum((positions - mean(positions, 2)).^2, 1)));
    ts.meanInputAmplitude(t) = mean(inputMetrics.amplitudeSum);
    ts.meanVisibleNeighbors(t) = mean(inputMetrics.visibleCount);
    ts.meanActiveNeighbors(t) = mean(inputMetrics.activeCount);
    ts.meanAngularCoverage(t) = mean(inputMetrics.angularCoverage);
    ts.meanInputConcentration(t) = mean(inputMetrics.inputConcentration);
    ts.meanFieldPeak(t) = mean(inputMetrics.fieldPeak);
    bumpDenominator = sum(positiveOutput, 1) + eps;
    ts.meanBumpConcentration(t) = mean(hypot(cx, cy)./bumpDenominator);
    ts.meanSpeed(t) = mean(speed);
    ts.meanInputMass(t) = mean(inputMetrics.inputMass);
    ts.perReceiverInputMass(t, :) = inputMetrics.inputMass(:).';
    if isfield(inputMetrics, 'hiddenWeightSum') && ~isempty(inputMetrics.hiddenWeightSum)
        ts.perReceiverHiddenWeight(t, :) = inputMetrics.hiddenWeightSum(:).';
    else
        ts.perReceiverHiddenWeight(t, :) = zeros(1, n);
    end
    ts.meanDistanceWeight(t) = mean(inputMetrics.distanceWeight);
    if cfg.scheduledOcclusion
        ts.meanPredictionError(t) = inputMetrics.meanPredictionAbsError;
        ts.meanPredictionSquaredError(t) = inputMetrics.meanPredictionSquaredError;
        ts.meanPredictionVariance(t) = inputMetrics.meanPredictionVariance;
    end
    if cfg.storeInputDiagnostics
        % Input-system diagnostics computed from the actually injected field (no
        % controller-visible extra information): the discrete circular mass and
        % the first circular harmonic of every receiver.
        discreteMass = sum(externalInput, 1)*2*pi/nS;
        harmonic = sum(externalInput .* exp(1i*alpha.'), 1);
        nonempty = discreteMass > 0;
        ts.meanInputMassExact(t) = mean(discreteMass);
        % Raw first circular moment divided by mass squared. NOTE: this is NOT
        % normalised, so its scale depends on the kernel convention (the mass-exact
        % wrapped kernel sums to h0; the historical truncated kernel sums to h0
        % times the grid integral, about 16). It must never be compared across
        % kernel modes. meanInputDirectionalInfo below is the scale-free companion.
        ts.meanInputHarmonic1(t) = mean(abs(harmonic(nonempty))./discreteMass(nonempty).^2);
        % Scale-free direction information: first circular harmonic divided by the
        % RAW grid sum of the field. Measured definition check (task3ProbeDirInfo):
        % a unit-mass wrapped normal with sigma=0.4 has sum=1 and |H1|=0.923, so
        % |H1|/sum = 0.923 and |H1|/(sum*2*pi/nS) = 14.7. Dividing by the
        % grid-integrated mass therefore inflates the value about 15-fold and can
        % exceed 1, which is impossible for a normalised kernel. This uses the raw
        % sum, giving a value in [0,1]: 1 for a point mass, falling toward 0 as the
        % covariance convolution flattens the kernel -- exactly the "lose
        % directional resolving power while keeping total mass" effect that
        % Section 7.2 requires be reported separately from mass.
        rawSum = sum(externalInput, 1);
        usable = rawSum > 0;
        ts.meanInputDirectionalInfo(t) = mean(abs(harmonic(usable))./rawSum(usable));
        if cfg.scheduledOcclusion && isfield(inputMetrics, 'inputMassRelError')
            % The plan requires the FFT-convolution mass error to be reported:
            % the clip at zero after the heat convolution is the only step that
            % can remove discrete mass, and hiding that would let a loss of
            % directional content be mistaken for a loss of total input.
            ts.meanInputMassRelError(t) = mean(inputMetrics.inputMassRelError);
        end
    end

    if t > cfg.burnIn
        adjacencySum = adjacencySum + adjacency;
        evaluationCount = evaluationCount + 1;
    end

    if cfg.storeTrajectory && nextRecord <= nRecords && t == recordSteps(nextRecord)
        xRecord(:, nextRecord) = positions(1, :).';
        yRecord(:, nextRecord) = positions(2, :).';
        headingRecord(:, nextRecord) = headings;
        nextRecord = nextRecord + 1;
    end
end

evaluationIndices = (cfg.burnIn + 1):T;
summary = struct();
metricNames = fieldnames(ts);
for idx = 1:numel(metricNames)
    name = metricNames{idx};
    val = ts.(name);
    % Only summarize 1D series (T x 1 vectors); exclude 2D receiver matrix diagnostics
    if ismatrix(val) && size(val, 2) == 1
        summary.(name) = mean(val(evaluationIndices), 'omitnan');
    end
end

if evaluationCount > 0
    meanAdjacency = adjacencySum/evaluationCount;
else
    meanAdjacency = zeros(n, n);
end
[summary.socialUnionLccFraction, summary.socialAlgebraicConnectivity] = ...
    graphSummary(meanAdjacency);

pilotState.finalNoiseStreamState = pilotState.noiseStream.State;
pilotState = rmfield(pilotState,'noiseStream');
result.pilot = pilotState;
result.config = cfg;
if cfg.scheduledOcclusion && cfg.storeObservationRecords
    result.observationRecords=observationRecords;
end
result.summary = summary;
if cfg.storeTimeSeries
    result.timeSeries = ts;
else
    result.timeSeries = [];
end
if cfg.storeTrajectory
    result.trajectory.time = recordSteps*cfg.dt;
    result.trajectory.x = xRecord;
    result.trajectory.y = yRecord;
    result.trajectory.headings = headingRecord;
else
    result.trajectory = [];
end
result.recordSteps = recordSteps;
% Post-initialization state at t=0, before any step was taken. Task 3 needs it
% to define the turn-cue base direction from the group heading at the start of
% the common warmup (Section 5.1).
result.initialPosition = initialPosition;
result.initialHeadings = initialHeadings;
result.initialNeuralStateUsed = initialNeuralState;
result.initialCueBaseAngle = angleOfCircularMeanLocal(initialHeadings);
result.cueTrace = cueTrace;
result.inputSamples=inputSamples;
result.lastSensorHeadings = sensorHeadings;
if cfg.scheduledOcclusion
    result.observationState = observationState;
end
result.finalPosition = positions;
result.finalHeadings = headings;
result.finalNeuralState = u;
result.finalPreferredDirections = alpha;
if cfg.storeMeanAdjacency
    result.meanAdjacency = meanAdjacency;
end
end

function [externalInput, adjacency, metrics] = buildSocialInput(positions, headings, alpha, cfg)
n = cfg.NAgents;
nS = cfg.NNeurons;
L = cfg.arenaSize;
kRequested = min(cfg.neighborCount, n - 1);

externalInput = zeros(nS, n);
adjacency = zeros(n, n);
metrics.visibleCount = zeros(n, 1);
metrics.activeCount = zeros(n, 1);
metrics.amplitudeSum = zeros(n, 1);
metrics.angularCoverage = zeros(n, 1);
metrics.inputConcentration = zeros(n, 1);
metrics.fieldPeak = zeros(n, 1);
metrics.inputMass = zeros(n,1);
metrics.distanceWeight = zeros(n,1);

baseAmplitude = cfg.totalSocialAttraction/n;
referenceTotal = (n - 1)*baseAmplitude;

for focal = 1:n
    candidates = [1:(focal - 1), (focal + 1):n];
    dx = positions(1, candidates) - positions(1, focal);
    dy = positions(2, candidates) - positions(2, focal);
    if ~(isfield(cfg, 'openArena') && cfg.openArena)
        dx = dx - L*round(dx/L);
        dy = dy - L*round(dy/L);
    end
    distances = hypot(dx, dy);
    bearings = mod(atan2(dy, dx), 2*pi);

    visible = visibilityMask(focal, candidates, positions, headings, ...
        dx, dy, distances, bearings, cfg);
    visibleLocal = find(visible);
    visibleCount = numel(visibleLocal);

    if strcmp(cfg.selectionPolicy, 'all')
        selectedLocal = visibleLocal;
    else
        selectedLocal = selectNeighbors( ...
            bearings(visible), distances(visible), ...
            min(kRequested, visibleCount), cfg.selectionPolicy);
        selectedLocal = visibleLocal(selectedLocal);
    end
    selectedAgents = candidates(selectedLocal);
    selectedBearings = bearings(selectedLocal);
    kActual = numel(selectedAgents);

    if strcmp(cfg.strengthMode, 'original')
        amplitude = baseAmplitude;
    elseif kActual > 0
        amplitude = referenceTotal/kActual;
    else
        amplitude = 0;
    end

    if kActual > 0
        angleDistance = abs(alpha(focal, :).'-selectedBearings);
        angleDistance = min(angleDistance, 2*pi-angleDistance);
        kernels=exp(-0.5*(angleDistance.^2)/(cfg.receptiveFieldWidth^2));
        if strcmp(cfg.distanceWeightMode,'equal')
            externalInput(:, focal)=amplitude*sum(kernels,2);
            metrics.distanceWeight(focal)=1;
        else
            [weights,normalized]=distanceWeight(distances(selectedLocal),cfg);
            if strcmp(cfg.strengthMode,'original')
                edgeAmplitudes=baseAmplitude*weights;
            else
                edgeAmplitudes=referenceTotal*normalized;
            end
            externalInput(:,focal)=sum(kernels.*edgeAmplitudes,2);
            metrics.distanceWeight(focal)=mean(weights);
        end
        adjacency(focal, selectedAgents) = 1;
    end

    metrics.visibleCount(focal) = visibleCount;
    metrics.activeCount(focal) = kActual;
    metrics.amplitudeSum(focal) = amplitude*kActual;
    if kActual>0 && ~strcmp(cfg.distanceWeightMode,'equal')
        metrics.amplitudeSum(focal)=sum(edgeAmplitudes);
    end
    metrics.inputMass(focal)=sum(externalInput(:,focal))*2*pi/nS;
    metrics.angularCoverage(focal) = angularCoverage(selectedBearings);
    if kActual > 0
        metrics.inputConcentration(focal) = abs(mean(exp(1i*selectedBearings)));
    end
    metrics.fieldPeak(focal) = max(externalInput(:, focal));
end

function visible = visibilityMask(focal, candidates, positions, headings, ...
        dx, dy, distances, bearings, cfg)
%VISIBILITYMASK Generate fresh bearing observations without range sensing.

visible = true(1, numel(candidates));
if strcmp(cfg.visibilityMode, 'bearing-fov')
    relativeBearing = wrapAngle(bearings-headings(focal));
    visible = abs(relativeBearing) <= cfg.fieldOfView/2;
    if cfg.detectionProbability < 1
        visible = visible & rand(1, numel(candidates)) < cfg.detectionProbability;
    end
elseif ~strcmp(cfg.visibilityMode, 'ideal')
    error('visibilityMode must be ideal or bearing-fov.');
end

if cfg.occlusionEnabled
    if cfg.agentRadius <= 0
        error('agentRadius must be positive when occlusionEnabled is true.');
    end
    blocked = false(1, numel(candidates));
    for target = 1:numel(candidates)
        if ~visible(target)
            continue;
        end
        for blocker = 1:numel(candidates)
            if blocker == target || ~visible(blocker)
                continue;
            end
            if distances(blocker) >= distances(target)
                continue;
            end
            targetDistance = distances(target);
            along = (dx(blocker)*dx(target) + dy(blocker)*dy(target)) / targetDistance;
            cross = abs(dx(blocker)*dy(target) - dy(blocker)*dx(target)) / targetDistance;
            if along > 0 && along < targetDistance && cross <= cfg.agentRadius
                blocked(target) = true;
                break;
            end
        end
    end
    visible(blocked) = false;
end
end

function angle = wrapAngle(angle)
angle = mod(angle + pi, 2*pi) - pi;
end
end

function selected = selectNeighbors(bearings, distances, k, policy)
nCandidates = numel(bearings);
if k <= 0
    selected = zeros(1, 0);
    return;
end
if k >= nCandidates
    selected = 1:nCandidates;
    return;
end

switch policy
    case 'random'
        selected = randperm(nCandidates, k);
    case 'nearest'
        [~, order] = sort(distances, 'ascend');
        selected = order(1:k);
    case 'balanced'
        [~, first] = min(distances);
        selected = first;
        available = true(1, nCandidates);
        available(first) = false;
        while numel(selected) < k
            candidates = find(available);
            separation = inf(size(candidates));
            for idx = 1:numel(candidates)
                delta = abs(bearings(candidates(idx))-bearings(selected));
                delta = min(delta, 2*pi-delta);
                separation(idx) = min(delta);
            end
            bestSeparation = max(separation);
            tied = candidates(abs(separation-bestSeparation) < 1e-12);
            [~, nearestTie] = min(distances(tied));
            chosen = tied(nearestTie);
            selected(end + 1) = chosen; %#ok<AGROW>
            available(chosen) = false;
        end
    otherwise
        error('Unknown selection policy: %s', policy);
end
end

function coverage = angularCoverage(bearings)
if numel(bearings) < 2
    coverage = 0;
    return;
end
ordered = sort(bearings);
gaps = [diff(ordered), 2*pi-ordered(end)+ordered(1)];
coverage = 1-max(gaps)/(2*pi);
end

function [meanPair, meanNearest] = periodicDistanceMetrics(positions, L, isOpen)
if nargin < 3, isOpen = false; end
dx = positions(1, :).'-positions(1, :);
dy = positions(2, :).'-positions(2, :);
if ~isOpen
    dx = dx-L*round(dx/L);
    dy = dy-L*round(dy/L);
end
distance = hypot(dx, dy);
n = size(distance, 1);
upper = distance(triu(true(n), 1));
meanPair = mean(upper);
distance(1:n+1:end) = inf;
meanNearest = mean(min(distance, [], 2));
end

function localOrder = topologicalLocalOrder(positions, unitX, unitY, L, k, isOpen)
%TOPOLOGICALLOCALORDER Match the paper's normalized VOP with k=5.
if nargin < 6, isOpen = false; end
n = size(positions, 2);
kLocal = min(k + 1, n);
dx = positions(1, :).'-positions(1, :);
dy = positions(2, :).'-positions(2, :);
if ~isOpen
    dx = dx-L*round(dx/L);
    dy = dy-L*round(dy/L);
end
distance = hypot(dx, dy);
distance(1:n+1:end) = 0;
localValues = zeros(n, 1);
for focal = 1:n
    [~, order] = sort(distance(focal, :), 'ascend');
    neighborhood = order(1:kLocal);
    localValues(focal) = hypot(sum(unitX(neighborhood)), ...
        sum(unitY(neighborhood)))/kLocal;
end
localOrder = mean(localValues);
end

function [lccFraction, algebraicConnectivity] = graphSummary(meanAdjacency)
n = size(meanAdjacency, 1);
unionAdjacency = (meanAdjacency + meanAdjacency.') > 0;
visited = false(n, 1);
largest = 0;
for start = 1:n
    if visited(start)
        continue;
    end
    queue = start;
    visited(start) = true;
    componentSize = 0;
    while ~isempty(queue)
        node = queue(1);
        queue(1) = [];
        componentSize = componentSize + 1;
        neighbors = find(unionAdjacency(node, :) & ~visited.');
        visited(neighbors) = true;
        queue = [queue, neighbors]; %#ok<AGROW>
    end
    largest = max(largest, componentSize);
end
lccFraction = largest/n;

weighted = (meanAdjacency + meanAdjacency.')/2;
laplacian = diag(sum(weighted, 2))-weighted;
eigenvalues = sort(real(eig(laplacian)));
if n >= 2
    algebraicConnectivity = eigenvalues(2);
else
    algebraicConnectivity = 0;
end
end

function mu = angleOfCircularMeanLocal(h)
h = h(:);
mu = atan2(mean(sin(h)), mean(cos(h)));
end

function cfg = validateConfig(cfg)
required = fieldnames(defaultEffectiveInputConfig());
defaults = defaultEffectiveInputConfig();
for idx = 1:numel(required)
    name = required{idx};
    if ~isfield(cfg, name)
        cfg.(name) = defaults.(name);
    end
end
% Task 3 fields are purely additive: only names absent from the base config
% are supplied, so a config that does not mention them keeps the original
% behaviour exactly.
base = fieldnames(defaultEffectiveInputConfig());
task3 = fieldnames(task3DefaultConfig());
t3defaults = task3DefaultConfig();
for idx = 1:numel(task3)
    name = task3{idx};
    if ~isfield(cfg, name) && ~any(strcmp(name, base))
        cfg.(name) = t3defaults.(name);
    end
end

validateattributes(cfg.NAgents, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(cfg.NNeurons, {'numeric'}, {'scalar','integer','>=',3});
validateattributes(cfg.NSteps, {'numeric'}, {'scalar','integer','>=',1});
validateattributes(cfg.burnIn, {'numeric'}, ...
    {'scalar','integer','>=',0,'<',cfg.NSteps});
validateattributes(cfg.neighborCount, {'numeric'}, ...
    {'scalar','integer','>=',0,'<=',cfg.NAgents - 1});
validateattributes(cfg.initialXY, {'numeric'}, {'2d'});
validateattributes(cfg.fieldOfView, {'numeric'}, {'scalar','>',0,'<=',2*pi});
validateattributes(cfg.detectionProbability, {'numeric'}, {'scalar','>=',0,'<=',1});
validateattributes(cfg.agentRadius, {'numeric'}, {'scalar','nonnegative'});
assert(islogical(cfg.turnExperiment) || (isnumeric(cfg.turnExperiment) && isscalar(cfg.turnExperiment)));
validateattributes(cfg.turnCueStartStep,{'numeric'},{'scalar','integer','>=',1});
validateattributes(cfg.turnCueDurationSteps,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.turnCueRate,{'numeric'},{'scalar','finite','real'});
validateattributes(cfg.turnCueAmplitude,{'numeric'},{'scalar','finite','nonnegative'});
if cfg.turnExperiment
    validateattributes(cfg.turnLeaderIds,{'numeric'},{'vector','integer','>=',1,'<=',cfg.NAgents});
    validateattributes(cfg.turnCueBaseAngles,{'numeric'},{'vector','numel',cfg.NAgents,'finite'});
    assert(numel(unique(cfg.turnLeaderIds))==numel(cfg.turnLeaderIds));
end
if cfg.scheduledOcclusion
 assert(cfg.allocentric && strcmp(cfg.strengthMode,'original') && strcmp(cfg.visibilityMode,'ideal') && ~cfg.occlusionEnabled);
 assert(strcmp(cfg.selectionPolicy,'fixed-id'));
 assert(any(strcmp(cfg.predictionMode,{'none','hold','mean','kf-mean','kalman','kalman-gm', ...
     'kf-conf-weight','kf-mixture-strip','kf-var-gated','kf-gm-damped', ...
     'hold-conf-weight','kf-global-massmatch','kf-relative-confidence'}))); % task2-relative-confidence-v1
 % Task 3 needs occlusionStartStep>=1 because a branch with zero occlusion
 % length is exactly the fresh reference. The original >=2 bound is preserved
 % for any caller that still passes a later start step.
 validateattributes(cfg.occlusionStartStep,{'numeric'},{'scalar','integer','>=',1});
 validateattributes(cfg.occlusionDurationSteps,{'numeric'},{'scalar','integer','nonnegative'});
 validateattributes(cfg.occlusionHiddenCount,{'numeric'},{'scalar','integer','>=',0,'<=',cfg.neighborCount});
 if ~isempty(cfg.occlusionHiddenIds)
     validateattributes(cfg.occlusionHiddenIds,{'numeric'},{'2d','integer','>=',1,'<=',cfg.NAgents});
     assert(size(cfg.occlusionHiddenIds,1)==cfg.NAgents);
     assert(size(cfg.occlusionHiddenIds,2)>=cfg.occlusionHiddenCount);
 end
 validateattributes(cfg.bearingProcessNoise,{'numeric'},{'scalar','finite','positive'});
 validateattributes(cfg.bearingGaussMarkovProcessNoise,{'numeric'},{'scalar','finite','positive'});
 validateattributes(cfg.bearingVelocityTau,{'numeric'},{'scalar','finite','positive'});
 validateattributes(cfg.bearingObservationVariance,{'numeric'},{'scalar','finite','positive'});
 validateattributes(cfg.bearingMeasurementStd,{'numeric'},{'scalar','finite','nonnegative'});
 validateattributes(cfg.predictionMaxAge,{'numeric'},{'scalar','finite','positive'});
 validateattributes(cfg.predictionMaxVariance,{'numeric'},{'scalar','finite','positive'});
end
validateattributes(cfg.absoluteStepOffset,{'numeric'},{'scalar','integer','nonnegative'});
assert(islogical(cfg.storeObservationRecords) && islogical(cfg.storeTrajectory) && ...
    islogical(cfg.storeTimeSeries) && islogical(cfg.storeInputDiagnostics));
if ~isempty(cfg.decompositionSteps)
    validateattributes(cfg.decompositionSteps,{'numeric'},{'vector','integer','>=',1,'<=',cfg.NSteps});
end
assert(any(strcmp(cfg.distanceWeightMode,{'equal','exponential'})));
validateattributes(cfg.distanceWeightScaleFraction,{'numeric'},{'scalar','real','finite','positive'});
if strcmp(cfg.distanceWeightMode,'exponential') && ...
        (~strcmp(cfg.visibilityMode,'ideal') || cfg.occlusionEnabled)
    error('Distance weights currently require ideal unobstructed range-bearing observations.');
end

policies = {'all','random','nearest','balanced'};
if cfg.scheduledOcclusion, policies{end+1}='fixed-id'; end
if ~any(strcmp(cfg.selectionPolicy, policies))
    error('selectionPolicy must be all, random, nearest, or balanced.');
end
modes = {'original','fixed-total'};
if ~any(strcmp(cfg.strengthMode, modes))
    error('strengthMode must be original or fixed-total.');
end
visibilityModes = {'ideal','bearing-fov'};
if ~any(strcmp(cfg.visibilityMode, visibilityModes))
    error('visibilityMode must be ideal or bearing-fov.');
end
if cfg.occlusionEnabled && ~strcmp(cfg.visibilityMode, 'bearing-fov')
    error('Occlusion requires visibilityMode=bearing-fov.');
end
if strcmp(cfg.selectionPolicy, 'all')
    cfg.neighborCount = cfg.NAgents - 1;
end

if ~isempty(cfg.initialXY) && ~isequal(size(cfg.initialXY), [2, cfg.NAgents])
    error('initialXY must be 2-by-NAgents.');
end
if ~isempty(cfg.initialHeadings) && numel(cfg.initialHeadings) ~= cfg.NAgents
    error('initialHeadings must contain NAgents values.');
end
if ~isempty(cfg.initialNeuralState) && ...
        ~isequal(size(cfg.initialNeuralState), [cfg.NNeurons, cfg.NAgents])
    error('initialNeuralState must be NNeurons-by-NAgents.');
end
if ~isempty(cfg.initialPreferredDirections) && ...
        ~isequal(size(cfg.initialPreferredDirections), [cfg.NAgents, cfg.NNeurons])
    error('initialPreferredDirections must be NAgents-by-NNeurons.');
end
end


