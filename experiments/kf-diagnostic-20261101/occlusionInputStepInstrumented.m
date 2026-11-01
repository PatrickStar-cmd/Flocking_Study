function [input,state,diag] = occlusionInputStepInstrumented(observation,alpha,state,opt)
%OCCLUSIONINPUTSTEPINSTRUMENTED Mechanically identical to stock occlusionInputStep,
% with added PASSIVE diagnostic recording (raw P, age, x, innovation, NIS, cutoffReason).
% Numerical paths, random state, and dynamics are 100% UNCHANGED.
%
% Cutoff reason codes:
%   -1: edge not selected
%   10: active fresh observation
%   11: active hold memory
%   12: active kf-mean memory
%   13: active kalman memory
%   20: truncated by 'none' mode during occlusion
%   21: truncated by age (> maxAge)
%   22: truncated by variance (P(1,1) > maxVariance)
%   23: truncated jointly by age and variance

n = size(observation.fresh, 1);
ns = size(alpha, 2);
if isempty(state)
    state.tracks = cell(n);
    state.bearing = nan(n);
    state.weight = nan(n);
    state.age = inf(n);
end

input = zeros(ns, n);
% Stock diagnostic fields (retained identically for 100% compatibility):
diag.active = zeros(n);
diag.variance = nan(n);
diag.center = nan(n);

% Passive diagnostic instrumentation fields:
diag.rawVariance = nan(n);
diag.rawP22 = nan(n);
diag.rawP12 = nan(n);
diag.rawP21 = nan(n);
diag.rawX = nan(n, n, 2);
diag.rawAge = nan(n);
diag.rawInnovation = nan(n);
diag.rawNis = nan(n);
diag.cutoffReason = repmat(-1, n, n);

for i = 1:n
    for j = find(observation.selected(i, :))
        fresh = observation.fresh(i, j);
        if fresh
            y = observation.bearing(i, j);
            w = observation.weight(i, j);
            assert(isfinite(y) && isfinite(w));
            state.bearing(i, j) = y;
            state.weight(i, j) = w;
            state.age(i, j) = 0;
        else
            assert(isnan(observation.bearing(i, j)) && isnan(observation.weight(i, j)));
            y = [];
            state.age(i, j) = state.age(i, j) + opt.dt;
        end

        if strcmp(opt.mode, 'kalman-gm')
            track = bearingTrackStepGaussMarkov(state.tracks{i, j}, y, opt.dt, opt.gmQ, opt.R, opt.velocityTau);
        else
            track = bearingTrackStep(state.tracks{i, j}, y, opt.dt, opt.q, opt.R);
        end
        state.tracks{i, j} = track;
        if isempty(track), continue; end

        % Record passive filter metrics regardless of whether edge is later active or truncated:
        diag.rawVariance(i, j) = track.P(1, 1);
        diag.rawP22(i, j) = track.P(2, 2);
        diag.rawP12(i, j) = track.P(1, 2);
        diag.rawP21(i, j) = track.P(2, 1);
        diag.rawX(i, j, 1) = track.x(1);
        diag.rawX(i, j, 2) = track.x(2);
        diag.rawAge(i, j) = state.age(i, j);
        diag.rawInnovation(i, j) = track.innovation;
        diag.rawNis(i, j) = track.nis;

        ageExceeded = (state.age(i, j) > opt.maxAge);
        varExceeded = (track.P(1, 1) > opt.maxVariance);

        if fresh
            diag.cutoffReason(i, j) = 10;
            center = y;
            variance = 0;
        elseif strcmp(opt.mode, 'none')
            diag.cutoffReason(i, j) = 20;
            continue;
        elseif state.age(i, j) > opt.maxAge
            if varExceeded && any(strcmp(opt.mode, {'mean','kf-mean','kalman','kalman-gm'}))
                diag.cutoffReason(i, j) = 23; % Joint age and variance cutoff
            else
                diag.cutoffReason(i, j) = 21; % Age cutoff
            end
            continue;
        elseif strcmp(opt.mode, 'hold')
            diag.cutoffReason(i, j) = 11;
            center = state.bearing(i, j);
            variance = 0;
        elseif strcmp(opt.mode, 'mean') || strcmp(opt.mode, 'kf-mean')
            if track.P(1, 1) > opt.maxVariance
                diag.cutoffReason(i, j) = 22; % Variance cutoff
                continue;
            end
            diag.cutoffReason(i, j) = 12;
            center = track.x(1);
            variance = 0;
        elseif strcmp(opt.mode, 'kalman') || strcmp(opt.mode, 'kalman-gm')
            if track.P(1, 1) > opt.maxVariance
                diag.cutoffReason(i, j) = 22; % Variance cutoff
                continue;
            end
            diag.cutoffReason(i, j) = 13;
            center = track.x(1);
            variance = track.P(1, 1);
        else
            error('Unknown prediction mode: %s', opt.mode);
        end

        diag.center(i, j) = center;
        delta = atan2(sin(alpha(i, :) - center), cos(alpha(i, :) - center));
        kernel = exp(-0.5 * (delta / opt.sigma).^2);
        if variance > 0
            % Circular heat convolution: preserves discrete mass, tends to uniform.
            frequency = [0:floor(ns/2), -ceil(ns/2)+1:-1];
            kernel = real(ifft(fft(kernel) .* exp(-0.5 * variance * frequency.^2)));
            kernel = max(0, kernel);
        end
        input(:, i) = input(:, i) + opt.h0 * state.weight(i, j) * kernel';
        diag.active(i, j) = 1;
        diag.variance(i, j) = track.P(1, 1);
    end
end
diag.age = state.age;
diag.weight = state.weight;
end
