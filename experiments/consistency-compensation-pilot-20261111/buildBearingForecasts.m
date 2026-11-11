function forecast = buildBearingForecasts(history, dt, maxSteps, params)
%BUILDBEARINGFORECASTS Predict missing neighbor bearings from pre-event data.
% history is H-by-N-by-K (radians), ending immediately before occlusion.
% Every returned field is maxSteps-by-N-by-K, wrapped to [0, 2*pi).
if nargin < 4 || isempty(params), params = struct(); end
assert(ndims(history) == 3 && size(history,1) >= 50, ...
    'history must be H-by-N-by-K with at least 50 samples.');
assert(isscalar(dt) && isfinite(dt) && dt > 0);
assert(isscalar(maxSteps) && maxSteps >= 1 && maxSteps == floor(maxSteps));
assert(all(isfinite(history(:))), 'Bearing history contains nonfinite values.');

alpha = option(params, 'alpha', 0.35);
beta = option(params, 'beta', 0.20);
phi = option(params, 'phi', 0.85);
q = option(params, 'q', 1e-5);
r = option(params, 'r', 1e-4);
historySteps = option(params, 'historySteps', 50);
assert(alpha > 0 && alpha <= 1 && beta > 0 && beta <= 1);
assert(phi > 0 && phi <= 1 && q > 0 && r > 0);
assert(historySteps >= 50 && historySteps == floor(historySteps));

H = size(history,1); N = size(history,2); K = size(history,3);
B = N*K; y = unwrap(reshape(history,H,B),[],1);
steps = (1:maxSteps)'*dt;
base = y(end,:);

forecast.Hold = shaped(repmat(base,maxSteps,1), maxSteps,N,K);

% A short least-squares line is the constant-angular-velocity baseline.
x = (-9:0)'*dt;
coef = [ones(10,1),x] \ y(end-9:end,:);
forecast.Linear = shaped(base + steps*coef(2,:),maxSteps,N,K);

% Holt's damped trend: the angular trend decays geometrically in blackout.
recent = y(end-min(H,historySteps)+1:end,:);
level = recent(1,:); trend = recent(2,:)-recent(1,:);
for t = 2:size(recent,1)
    oldLevel = level;
    level = alpha*recent(t,:) + (1-alpha)*(level+phi*trend);
    trend = beta*(level-oldLevel) + (1-beta)*phi*trend;
end
if phi == 1
    trendSum = (1:maxSteps)';
else
    trendSum = phi*(1-phi.^(1:maxSteps)')/(1-phi);
end
forecast.DampedHolt = shaped(level + trendSum*trend,maxSteps,N,K);

% Savitzky-Golay derivative at the last observed sample; extrapolate the
% tangent rather than an unbounded quadratic trajectory.
x = (-10:0)'*dt;
coef = [ones(11,1),x,x.^2] \ y(end-10:end,:);
forecast.SG = shaped(base + steps*coef(2,:),maxSteps,N,K);

% Circular constant-velocity Kalman filter with white angular acceleration.
% Unwrapped observations keep the innovation continuous across 2*pi.
xAngle = recent(1,:);
xRate = (recent(2,:)-recent(1,:))/dt;
P11 = r*ones(1,B); P12 = zeros(1,B); P22 = ones(1,B);
Q11 = q*dt^3/3; Q12 = q*dt^2/2; Q22 = q*dt;
for t = 2:size(recent,1)
    priorAngle = xAngle + dt*xRate;
    priorP11 = P11 + 2*dt*P12 + dt^2*P22 + Q11;
    priorP12 = P12 + dt*P22 + Q12;
    priorP22 = P22 + Q22;
    innovation = recent(t,:) - priorAngle;
    K1 = priorP11./(priorP11+r);
    K2 = priorP12./(priorP11+r);
    xAngle = priorAngle + K1.*innovation;
    xRate = xRate + K2.*innovation;
    P11 = (1-K1).*priorP11;
    P12 = (1-K1).*priorP12;
    P22 = priorP22 - K2.*priorP12;
end
forecast.KF = shaped(xAngle + steps*xRate,maxSteps,N,K);
end

function value = option(options, name, fallback)
if isfield(options,name), value = options.(name); else, value = fallback; end
end

function angles = shaped(values, steps, N, K)
angles = reshape(mod(values,2*pi),steps,N,K);
end
