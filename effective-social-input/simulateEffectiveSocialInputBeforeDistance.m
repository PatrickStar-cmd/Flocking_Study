function result = simulateEffectiveSocialInputBeforeDistance(cfg)
%SIMULATEEFFECTIVESOCIALINPUT Allocentric flocking with controlled bearing inputs.
%
% This implementation follows the neural-field update in Salahshour and
% Couzin (2025), while exposing the number and angular selection of social
% bearing inputs. It intentionally does not add range sensing or prediction.

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
alpha0 = linspace(0, 2*pi, nS + 1);
alpha0 = alpha0(1:end-1);
alpha = mod(repmat(alpha0, n, 1) + headings, 2*pi);

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

adjacencySum = zeros(n, n);
evaluationCount = 0;
cosAlpha = cos(alpha).';
sinAlpha = sin(alpha).';

for t = 1:T
    [externalInput, adjacency, inputMetrics] = buildSocialInput( ...
        positions, headings, alpha, cfg);

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
    positions(1, :) = mod(positions(1, :) + cfg.dt*velocityX, L);
    positions(2, :) = mod(positions(2, :) + cfg.dt*velocityY, L);

    speed = hypot(velocityX, velocityY);
    unitX = zeros(1, n);
    unitY = zeros(1, n);
    nonzero = speed > 1e-12;
    unitX(nonzero) = velocityX(nonzero)./speed(nonzero);
    unitY(nonzero) = velocityY(nonzero)./speed(nonzero);
    ts.globalOrder(t) = hypot(sum(unitX), sum(unitY))/n;
    ts.localOrder(t) = topologicalLocalOrder(positions, unitX, unitY, L, 5);

    [meanPair, meanNearest] = periodicDistanceMetrics(positions, L);
    ts.meanPairDistanceNorm(t) = meanPair/L;
    ts.meanNearestDistanceNorm(t) = meanNearest/L;
    ts.meanInputAmplitude(t) = mean(inputMetrics.amplitudeSum);
    ts.meanVisibleNeighbors(t) = mean(inputMetrics.visibleCount);
    ts.meanActiveNeighbors(t) = mean(inputMetrics.activeCount);
    ts.meanAngularCoverage(t) = mean(inputMetrics.angularCoverage);
    ts.meanInputConcentration(t) = mean(inputMetrics.inputConcentration);
    ts.meanFieldPeak(t) = mean(inputMetrics.fieldPeak);
    bumpDenominator = sum(positiveOutput, 1) + eps;
    ts.meanBumpConcentration(t) = mean(hypot(cx, cy)./bumpDenominator);
    ts.meanSpeed(t) = mean(speed);

    if t > cfg.burnIn
        adjacencySum = adjacencySum + adjacency;
        evaluationCount = evaluationCount + 1;
    end

    if nextRecord <= nRecords && t == recordSteps(nextRecord)
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
    summary.(name) = mean(ts.(name)(evaluationIndices), 'omitnan');
end

if evaluationCount > 0
    meanAdjacency = adjacencySum/evaluationCount;
else
    meanAdjacency = zeros(n, n);
end
[summary.socialUnionLccFraction, summary.socialAlgebraicConnectivity] = ...
    graphSummary(meanAdjacency);

result.config = cfg;
result.summary = summary;
result.timeSeries = ts;
result.trajectory.time = recordSteps*cfg.dt;
result.trajectory.x = xRecord;
result.trajectory.y = yRecord;
result.trajectory.headings = headingRecord;
result.finalPosition = positions;
result.finalHeadings = headings;
result.finalNeuralState = u;
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

baseAmplitude = cfg.totalSocialAttraction/n;
referenceTotal = (n - 1)*baseAmplitude;

for focal = 1:n
    candidates = [1:(focal - 1), (focal + 1):n];
    dx = positions(1, candidates) - positions(1, focal);
    dy = positions(2, candidates) - positions(2, focal);
    dx = dx - L*round(dx/L);
    dy = dy - L*round(dy/L);
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
        externalInput(:, focal) = amplitude*sum( ...
            exp(-0.5*(angleDistance.^2)/(cfg.receptiveFieldWidth^2)), 2);
        adjacency(focal, selectedAgents) = 1;
    end

    metrics.visibleCount(focal) = visibleCount;
    metrics.activeCount(focal) = kActual;
    metrics.amplitudeSum(focal) = amplitude*kActual;
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

function [meanPair, meanNearest] = periodicDistanceMetrics(positions, L)
dx = positions(1, :).'-positions(1, :);
dy = positions(2, :).'-positions(2, :);
dx = dx-L*round(dx/L);
dy = dy-L*round(dy/L);
distance = hypot(dx, dy);
n = size(distance, 1);
upper = distance(triu(true(n), 1));
meanPair = mean(upper);
distance(1:n+1:end) = inf;
meanNearest = mean(min(distance, [], 2));
end

function localOrder = topologicalLocalOrder(positions, unitX, unitY, L, k)
%TOPOLOGICALLOCALORDER Match the paper's normalized VOP with k=5.
n = size(positions, 2);
kLocal = min(k + 1, n);
dx = positions(1, :).'-positions(1, :);
dy = positions(2, :).'-positions(2, :);
dx = dx-L*round(dx/L);
dy = dy-L*round(dy/L);
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

function cfg = validateConfig(cfg)
required = fieldnames(defaultEffectiveInputConfig());
defaults = defaultEffectiveInputConfig();
for idx = 1:numel(required)
    name = required{idx};
    if ~isfield(cfg, name)
        cfg.(name) = defaults.(name);
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

policies = {'all','random','nearest','balanced'};
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
end
