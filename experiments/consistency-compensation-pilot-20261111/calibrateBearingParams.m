function params = calibrateBearingParams(histories, dt, horizons)
%CALIBRATEBEARINGPARAMS Tune two forecast models on pre-occlusion histories.
% histories is a cell array of H-by-N-by-K bearing histories from training
% seeds only. Targets are held-out future samples within each history.
if nargin < 3 || isempty(horizons), horizons = [30,100]; end
assert(iscell(histories) && ~isempty(histories));
assert(all(horizons >= 1 & horizons == floor(horizons)));

holtGrid = allcomb([0.2,0.5,0.8],[0.05,0.2,0.5],[0.6,0.85,0.97]);
kfGrid = allcomb([1e-6,1e-5,1e-4],[1e-6,1e-4,1e-3]);
bestHolt = inf; bestKF = inf;
params = struct('alpha',0.35,'beta',0.2,'phi',0.85, ...
    'q',1e-5,'r',1e-4,'historySteps',50);

for c = 1:size(holtGrid,1)
    candidate = params;
    candidate.alpha = holtGrid(c,1);
    candidate.beta = holtGrid(c,2);
    candidate.phi = holtGrid(c,3);
    mse = score(histories,dt,horizons,candidate,'DampedHolt');
    if mse < bestHolt
        bestHolt = mse;
        params.alpha = candidate.alpha;
        params.beta = candidate.beta;
        params.phi = candidate.phi;
    end
end
for c = 1:size(kfGrid,1)
    candidate = params;
    candidate.q = kfGrid(c,1);
    candidate.r = kfGrid(c,2);
    mse = score(histories,dt,horizons,candidate,'KF');
    if mse < bestKF
        bestKF = mse;
        params.q = candidate.q;
        params.r = candidate.r;
    end
end
params.calibrationHoltMSE = bestHolt;
params.calibrationKFMSE = bestKF;
params.calibrationSeeds = numel(histories);
params.calibrationHorizons = horizons;
end

function mse = score(histories,dt,horizons,params,method)
sumSquared = 0; count = 0;
for si = 1:numel(histories)
    y = histories{si}; H = size(y,1);
    assert(ndims(y) == 3 && H >= 150);
    anchors = unique([50,max(50,round(H*[0.35,0.50,0.65,0.80]))]);
    for a = anchors
        valid = horizons(a+horizons <= H);
        if isempty(valid), continue; end
        f = buildBearingForecasts(y(1:a,:,:),dt,max(valid),params);
        p = f.(method);
        for h = valid
            predicted = reshape(p(h,:,:),[],1);
            truth = reshape(y(a+h,:,:),[],1);
            error = atan2(sin(predicted-truth),cos(predicted-truth));
            sumSquared = sumSquared + sum(error.^2);
            count = count + numel(error);
        end
    end
end
assert(count > 0, 'No valid calibration anchors.');
mse = sumSquared/count;
end

function grid = allcomb(varargin)
v = cell(1,nargin);
[v{:}] = ndgrid(varargin{:});
grid = zeros(numel(v{1}),nargin);
for j = 1:nargin, grid(:,j) = v{j}(:); end
end
