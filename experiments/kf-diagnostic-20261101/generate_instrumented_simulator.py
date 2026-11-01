#!/usr/bin/env python3
"""
generate_instrumented_simulator.py

Generates passive instrumented simulator copies in experiments/kf-diagnostic-20261101/
from stock sources in effective-social-input/.

Rules:
1. No changes to numerical paths, random number calls, or core flocking dynamics.
2. Passive recording only: uRecord, positionsRecord, headingsRecord, stateFingerprints,
   preOcclusionCheckpoint, rawVariance, rawP22, rawX, rawAge, rawInnovation, rawNis, cutoffReason.
3. Completely avoids shadowing stock functions by using distinct filenames:
   - simulateEffectiveSocialInputInstrumented.m
   - scheduledOcclusionInputInstrumented.m
   - occlusionInputStepInstrumented.m
4. Computes SHA256 hashes and generates unified diff.
"""

import os
import hashlib
import difflib
import json

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC_DIR = os.path.join(BASE_DIR, "effective-social-input")
OUT_DIR = os.path.abspath(os.path.dirname(__file__))

def compute_sha256(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def make_diff(original_text, new_text, original_name, new_name):
    orig_lines = original_text.splitlines(keepends=True)
    new_lines = new_text.splitlines(keepends=True)
    diff = difflib.unified_diff(orig_lines, new_lines, fromfile=original_name, tofile=new_name)
    return "".join(diff)

def main():
    print(f"Base Directory: {BASE_DIR}")
    print(f"Source Directory: {SRC_DIR}")
    print(f"Output Directory: {OUT_DIR}")

    stock_files = {
        "simulator": os.path.join(SRC_DIR, "simulateEffectiveSocialInput.m"),
        "scheduled_sensor": os.path.join(SRC_DIR, "scheduledOcclusionInput.m"),
        "step_sensor": os.path.join(SRC_DIR, "occlusionInputStep.m"),
        "bearing_track": os.path.join(SRC_DIR, "bearingTrackStep.m"),
    }

    stock_hashes = {}
    for name, path in stock_files.items():
        if os.path.exists(path):
            stock_hashes[name] = {
                "path": os.path.relpath(path, BASE_DIR).replace("\\", "/"),
                "sha256": compute_sha256(path)
            }
        else:
            raise FileNotFoundError(f"Stock file missing: {path}")

    # Generate instrumented files
    out_step_path = os.path.join(OUT_DIR, "occlusionInputStepInstrumented.m")
    out_sensor_path = os.path.join(OUT_DIR, "scheduledOcclusionInputInstrumented.m")
    out_sim_path = os.path.join(OUT_DIR, "simulateEffectiveSocialInputInstrumented.m")

    # Read stock contents
    with open(stock_files["step_sensor"], "r", encoding="utf-8") as f:
        stock_step_code = f.read()
    with open(stock_files["scheduled_sensor"], "r", encoding="utf-8") as f:
        stock_sensor_code = f.read()
    with open(stock_files["simulator"], "r", encoding="utf-8") as f:
        stock_sim_code = f.read()

    # Write instrumented occlusionInputStepInstrumented.m
    inst_step_code = '''function [input,state,diag] = occlusionInputStepInstrumented(observation,alpha,state,opt)
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
'''

    # Write instrumented scheduledOcclusionInputInstrumented.m
    inst_sensor_code = '''function [input,state,adj,m,record] = scheduledOcclusionInputInstrumented(position,alpha,state,cfg,step)
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
'''

    # Write instrumented simulateEffectiveSocialInputInstrumented.m
    # We mechanically replicate simulateEffectiveSocialInput.m but:
    # 1. Rename to simulateEffectiveSocialInputInstrumented
    # 2. Call scheduledOcclusionInputInstrumented instead of scheduledOcclusionInput
    # 3. Add uRecord, positionsRecord, headingsRecord, stateFingerprints, preOcclusionCheckpoint
    # 4. Save them in result.
    inst_sim_code = stock_sim_code.replace(
        "function result = simulateEffectiveSocialInput(cfg)",
        "function result = simulateEffectiveSocialInputInstrumented(cfg)\n"
        "%SIMULATEEFFECTIVESOCIALINPUTINSTRUMENTED Passive instrumented version of simulateEffectiveSocialInput.\n"
        "% Retains exact identical dynamics while passively recording full uRecord, positionsRecord,\n"
        "% headingsRecord, stateFingerprints, preOcclusionCheckpoint, and full observation records."
    )
    inst_sim_code = inst_sim_code.replace(
        "[externalInput,observationState,adjacency,inputMetrics,observationRecords{t}]=\n            scheduledOcclusionInput(positions,alpha,observationState,cfg,t);",
        "[externalInput,observationState,adjacency,inputMetrics,observationRecords{t}]= ...\n            scheduledOcclusionInputInstrumented(positions,alpha,observationState,cfg,t);"
    )
    inst_sim_code = inst_sim_code.replace(
        "[externalInput,observationState,adjacency,inputMetrics,observationRecords{t}]= ...\n            scheduledOcclusionInput(positions,alpha,observationState,cfg,t);",
        "[externalInput,observationState,adjacency,inputMetrics,observationRecords{t}]= ...\n            scheduledOcclusionInputInstrumented(positions,alpha,observationState,cfg,t);"
    )

    # Insert allocations before simulation loop
    allocation_anchor = "cosAlpha = cos(alpha).';\nsinAlpha = sin(alpha).';"
    allocation_insert = (
        "cosAlpha = cos(alpha).';\n"
        "sinAlpha = sin(alpha).';\n"
        "\n"
        "% --- PASSIVE DIAGNOSTIC ALLOCATION ---\n"
        "uRecord = zeros(nS, n, T);\n"
        "positionsRecord = zeros(2, n, T);\n"
        "headingsRecord = zeros(n, T);\n"
        "stateFingerprints = zeros(T, 1);\n"
        "preOcclusionCheckpoint = struct();\n"
        "hasRecordedCheckpoint = false;\n"
    )
    if allocation_anchor in inst_sim_code:
        inst_sim_code = inst_sim_code.replace(allocation_anchor, allocation_insert)
    else:
        raise RuntimeError("Allocation anchor not found in stock simulator code.")

    # Insert recording at each step in loop
    loop_anchor = "for t = 1:T\n"
    loop_insert = (
        "for t = 1:T\n"
        "    % --- PASSIVE DIAGNOSTIC STEP RECORDING (pre-dynamics of step t) ---\n"
        "    positionsRecord(:, :, t) = positions;\n"
        "    headingsRecord(:, t) = headings;\n"
        "    uRecord(:, :, t) = u;\n"
        "    % Deterministic fingerprint of (positions, headings, u)\n"
        "    stateVec = [positions(:); headings(:); u(:)];\n"
        "    stateFingerprints(t) = sum(stateVec .* sin((1:numel(stateVec))' * 0.1));\n"
        "\n"
        "    if cfg.scheduledOcclusion && t == cfg.occlusionStartStep && ~hasRecordedCheckpoint\n"
        "        preOcclusionCheckpoint.step = t;\n"
        "        preOcclusionCheckpoint.stateAfterStep = t - 1;\n"
        "        preOcclusionCheckpoint.rngState = rng;\n"
        "        preOcclusionCheckpoint.positions = positions;\n"
        "        preOcclusionCheckpoint.headings = headings;\n"
        "        preOcclusionCheckpoint.u = u;\n"
        "        preOcclusionCheckpoint.alpha = alpha;\n"
        "        if exist('observationState', 'var') && ~isempty(observationState)\n"
        "            preOcclusionCheckpoint.observationState = observationState;\n"
        "        end\n"
        "        hasRecordedCheckpoint = true;\n"
        "    end\n"
    )
    if loop_anchor in inst_sim_code:
        inst_sim_code = inst_sim_code.replace(loop_anchor, loop_insert, 1)
    else:
        raise RuntimeError("Loop anchor not found in stock simulator code.")

    # Insert fields into result struct
    result_anchor = "result.finalPreferredDirections = alpha;"
    result_insert = (
        "result.finalPreferredDirections = alpha;\n"
        "% --- PASSIVE DIAGNOSTIC RESULTS ---\n"
        "result.uRecord = uRecord;\n"
        "result.positionsRecord = positionsRecord;\n"
        "result.headingsRecord = headingsRecord;\n"
        "result.stateFingerprints = stateFingerprints;\n"
        "result.preOcclusionCheckpoint = preOcclusionCheckpoint;\n"
        "if cfg.scheduledOcclusion && exist('observationState', 'var') && isfield(observationState, 'sensorNoise')\n"
        "    result.sensorNoise = observationState.sensorNoise;\n"
        "end"
    )
    if result_anchor in inst_sim_code:
        inst_sim_code = inst_sim_code.replace(result_anchor, result_insert)
    else:
        raise RuntimeError("Result anchor not found in stock simulator code.")

    # Write out files
    with open(out_step_path, "w", encoding="utf-8") as f:
        f.write(inst_step_code)
    with open(out_sensor_path, "w", encoding="utf-8") as f:
        f.write(inst_sensor_code)
    with open(out_sim_path, "w", encoding="utf-8") as f:
        f.write(inst_sim_code)

    print(f"Generated: {out_step_path}")
    print(f"Generated: {out_sensor_path}")
    print(f"Generated: {out_sim_path}")

    # Generate unified diffs
    diff_step = make_diff(stock_step_code, inst_step_code, "effective-social-input/occlusionInputStep.m", "experiments/kf-diagnostic-20261101/occlusionInputStepInstrumented.m")
    diff_sensor = make_diff(stock_sensor_code, inst_sensor_code, "effective-social-input/scheduledOcclusionInput.m", "experiments/kf-diagnostic-20261101/scheduledOcclusionInputInstrumented.m")
    diff_sim = make_diff(stock_sim_code, inst_sim_code, "effective-social-input/simulateEffectiveSocialInput.m", "experiments/kf-diagnostic-20261101/simulateEffectiveSocialInputInstrumented.m")

    full_diff = f"{diff_step}\n\n{diff_sensor}\n\n{diff_sim}\n"
    diff_path = os.path.join(OUT_DIR, "instrumentation.diff")
    with open(diff_path, "w", encoding="utf-8") as f:
        f.write(full_diff)
    print(f"Generated Diff: {diff_path}")

    # Record hashes
    manifest = {
        "stock_sources": stock_hashes,
        "instrumented_copies": {
            "occlusionInputStepInstrumented": {
                "path": "experiments/kf-diagnostic-20261101/occlusionInputStepInstrumented.m",
                "sha256": compute_sha256(out_step_path)
            },
            "scheduledOcclusionInputInstrumented": {
                "path": "experiments/kf-diagnostic-20261101/scheduledOcclusionInputInstrumented.m",
                "sha256": compute_sha256(out_sensor_path)
            },
            "simulateEffectiveSocialInputInstrumented": {
                "path": "experiments/kf-diagnostic-20261101/simulateEffectiveSocialInputInstrumented.m",
                "sha256": compute_sha256(out_sim_path)
            }
        },
        "diff": {
            "path": "experiments/kf-diagnostic-20261101/instrumentation.diff",
            "sha256": compute_sha256(diff_path)
        }
    }

    hash_path = os.path.join(OUT_DIR, "source_hashes.json")
    with open(hash_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Saved Hashes: {hash_path}")

if __name__ == "__main__":
    main()
