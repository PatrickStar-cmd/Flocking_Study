function runConsistencyCompensation2(mode,seedOverride)
% Paired intermittent-occlusion consistency experiment, with a locked gate.
% Modes: plan, smoke, development, calibrate, test.
if nargin < 1, mode = 'plan'; end
if nargin < 2, seedOverride = []; end
root = fileparts(mfilename('fullpath'));
cfg = jsondecode(fileread(fullfile(root, 'config.json')));
validateConfig(cfg);
seeds = cfg.seed_first + (0:cfg.seed_count-1);
if strcmp(mode, 'plan')
    fprintf('%d seeds: %d development, %d held-out test.\n', ...
        numel(seeds), cfg.development_seeds, numel(seeds)-cfg.development_seeds);
    for j = 1:numel(cfg.schedule_candidates)
        s = cfg.schedule_candidates(j);
        fprintf('%s: %d blackout + %d refresh, %d cycles.\n', ...
            s.name, s.occlusion_steps, s.refresh_steps, s.cycles);
    end
    return;
end
assert(any(strcmp(mode, {'smoke','development','calibrate','test'})), 'Unknown mode.');
if strcmp(mode, 'smoke')
    runSmoke(seeds(1), cfg, root);
    return;
end
if strcmp(mode, 'development')
    devSeeds = seeds(1:cfg.development_seeds);
    if ~isempty(seedOverride)
        assert(all(ismember(seedOverride,devSeeds)) && ...
            numel(unique(seedOverride)) == numel(seedOverride), ...
            'Development override must contain unique development seeds.');
        devSeeds = seedOverride(:).';
    end
    started = tic;
    for j = 1:numel(cfg.schedule_candidates)
        s = cfg.schedule_candidates(j);
        out = fullfile(root, 'results', 'development', s.name);
        freezeConfig(root, out);
        for si = 1:numel(devSeeds)
            seed = devSeeds(si); folder = fullfile(out, sprintf('seed-%d', seed));
            if isComplete(folder, seed, s.name, cfg.remaining_counts, {'Fresh','Drop'})
                fprintf('Development reused %s seed %d.\n', s.name, seed);
                continue;
            end
            runSeed(seed, cfg, s, {'Fresh','Drop'}, struct(), out, 'development');
            fprintf('Development %s seed %d (%d/%d), %.1f min.\n', ...
                s.name, seed, si, numel(devSeeds), toc(started)/60);
        end
        if isempty(seedOverride), writeCompletion(out, cfg.development_seeds, s); end
    end
    return;
end
gateFile = fullfile(root, 'results', 'development', 'analysis', 'gate.json');
assert(exist(gateFile, 'file') == 2, 'Run Python development analysis before calibrate/test.');
gate = jsondecode(fileread(gateFile));
assert(isfield(gate,'status') && strcmpi(gate.status,'pass'), ...
    'Development gate did not pass; test run is blocked.');
assert(isfield(gate,'selected_schedule') && ischar(gate.selected_schedule), ...
    'Gate has no selected schedule.');
s = findSchedule(cfg, gate.selected_schedule);
paramFile = fullfile(root, 'results', 'test', 'forecast_params.json');
if strcmp(mode, 'calibrate') || ~exist(paramFile, 'file')
    histories = cell(cfg.development_seeds,1);
    for si = 1:cfg.development_seeds
        [~,histories{si}] = makeCheckpoint(seeds(si),cfg);
    end
    params = calibrateBearingParams(histories, cfg.model.dt, ...
        unique([round(s.occlusion_steps/2),s.occlusion_steps]));
    params.schedule = s.name;
    out = fullfile(root, 'results', 'test');
    freezeConfig(root,out);
    writeJson(paramFile,params);
    fprintf('Calibrated parameters for %s on %d development seeds.\n', ...
        s.name, cfg.development_seeds);
    if strcmp(mode,'calibrate'), return; end
end
params = jsondecode(fileread(paramFile));
assert(strcmp(params.schedule,s.name), 'Stored forecast parameters have another schedule.');
out = fullfile(root, 'results', 'test');
freezeConfig(root,out);
methods = cfg.methods;
if isstring(methods), methods = cellstr(methods); end
testSeeds = seeds(cfg.development_seeds+1:end);
if ~isempty(seedOverride)
    assert(all(ismember(seedOverride,testSeeds)) && ...
        numel(unique(seedOverride)) == numel(seedOverride), ...
        'Test override must contain unique held-out seeds.');
    testSeeds = seedOverride(:).';
end
started = tic;
for si = 1:numel(testSeeds)
    seed = testSeeds(si); folder = fullfile(out,sprintf('seed-%d',seed));
    if isComplete(folder,seed,s.name,cfg.remaining_counts,methods)
        fprintf('Test reused seed %d.\n',seed);
        continue;
    end
    runSeed(seed,cfg,s,methods,params,out,'test');
    fprintf('Test seed %d (%d/%d), %.1f min.\n',seed, ...
        si,numel(testSeeds),toc(started)/60);
end
if isempty(seedOverride), writeCompletion(out,numel(testSeeds),s); end
end

function validateConfig(cfg)
assert(cfg.model.N == 20 && cfg.model.full_neighbors == cfg.model.N-1);
assert(cfg.seed_count == 30 && cfg.development_seeds == 6);
assert(cfg.history_steps == 50 && cfg.pre_steps >= cfg.history_steps);
assert(cfg.kappa == 1 && cfg.observation_noise_sd == 0);
assert(isequal(sort(cfg.remaining_counts(:)).',[5,8]));
assert(cfg.primary_remaining == 5 && cfg.tail_steps == 1000);
assert(numel(cfg.schedule_candidates) == 2);
for j = 1:numel(cfg.schedule_candidates)
    s = cfg.schedule_candidates(j);
    assert(s.refresh_steps == cfg.history_steps);
    assert(s.cycles*(s.occlusion_steps+s.refresh_steps) == 6000);
end
end

function s = findSchedule(cfg,name)
matches = strcmp({cfg.schedule_candidates.name},name);
assert(sum(matches) == 1, 'Gate selected a schedule absent from config.');
s = cfg.schedule_candidates(matches);
end

function freezeConfig(root,out)
if ~exist(out,'dir'), mkdir(out); end
source = fullfile(root,'config.json'); target = fullfile(out,'config.json');
if exist(target,'file')
    assert(strcmp(fileread(source),fileread(target)), 'Frozen result config differs.');
else
    copyfile(source,target);
end
end

function writeJson(path,obj)
fid = fopen(path,'w'); assert(fid > 0,'Cannot write %s.',path);
fprintf(fid,'%s\n',jsonencode(obj,'PrettyPrint',true)); fclose(fid);
end

function writeCompletion(out,nSeeds,s)
obj = struct('status','complete','seed_count',nSeeds,'schedule',s.name, ...
    'post_steps',s.cycles*(s.occlusion_steps+s.refresh_steps), ...
    'finished',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
writeJson(fullfile(out,'completion.json'),obj);
end

function tf = isComplete(folder,seed,schedule,remainingCounts,methods)
f = fullfile(folder,'summary.csv'); g = fullfile(folder,'traces.mat');
tf = exist(f,'file') == 2 && exist(g,'file') == 2;
if ~tf, return; end
try
    t = readtable(f);
    expected = numel(remainingCounts)*numel(methods);
    tf = height(t) == expected && all(t.seed == seed) && ...
        all(strcmp(t.schedule,schedule));
    if ~tf, return; end
    for r = remainingCounts(:).'
        for m = methods(:).'
            tf = tf && sum(t.remaining == r & strcmp(t.method,m{1})) == 1;
        end
    end
catch
    tf = false;
end
end

function runSeed(seed,cfg,s,methods,params,out,split)
[checkpoint,history,c,order] = makeCheckpoint(seed,cfg);
recent = history(end-cfg.history_steps+1:end,:,:);
rows = struct([]); traces = struct([]);
fresh = simulateBranch(checkpoint,recent,c,true(c.N,c.full_neighbors), ...
    'Fresh',params,s,true);
for r = cfg.remaining_counts(:).'
    visible = visibleMask(order,r,c);
    S0 = initialSocialInput(checkpoint,visible,c,cfg.kappa);
    assert(S0 >= 0 && S0 <= r*c.single_edge_weight+1e-12);
    for mi = 1:numel(methods)
        method = methods{mi};
        if strcmp(method,'Fresh')
            result = fresh;
        else
            result = simulateBranch(checkpoint,recent,c,visible,method,params,s,true);
        end
        row = summarizeBranch(seed,split,s.name,r,method,S0,result.series,cfg.tail_steps);
        rows = appendStruct(rows,row);
        traces = appendStruct(traces,struct('remaining',r,'method',method,'series',result.series));
    end
end
folder = fullfile(out,sprintf('seed-%d',seed));
if ~exist(folder,'dir'), mkdir(folder); end
writetable(struct2table(rows),fullfile(folder,'summary.tmp.csv'));
save(fullfile(folder,'traces.tmp.mat'),'traces','checkpoint','order','params','-v7');
movefile(fullfile(folder,'traces.tmp.mat'),fullfile(folder,'traces.mat'),'f');
movefile(fullfile(folder,'summary.tmp.csv'),fullfile(folder,'summary.csv'),'f');
end

function row = summarizeBranch(seed,split,schedule,r,method,S0,series,tailSteps)
total = size(series,1); tail = false(total,1);
tail(total-min(tailSteps,total)+1:total) = true;
occluded = series(:,10) == 1; occTail = occluded & tail;
assert(any(occTail) && any(occluded));
row = struct('seed',seed,'split',split,'schedule',schedule, ...
    'remaining',r,'method',method,'S0',S0, ...
    'C_occ_tail',mean(series(occTail,3)), ...
    'GO_occ_tail',mean(series(occTail,1)), ...
    'LO_occ_tail',mean(series(occTail,2)), ...
    'C_all_tail',mean(series(tail,3)), ...
    'GO_all_tail',mean(series(tail,1)), ...
    'LO_all_tail',mean(series(tail,2)), ...
    'C_occ_all',mean(series(occluded,3)), ...
    'GO_occ_all',mean(series(occluded,1)), ...
    'LO_occ_all',mean(series(occluded,2)), ...
    'bearingMAE',mean(series(occluded,8),'omitnan'), ...
    'inputRMSE',sqrt(mean(series(occluded,9).^2)));
end

function x = appendStruct(x,row)
if isempty(x), x = row; else, x(end+1) = row; end
end

function [checkpoint,history,c,order] = makeCheckpoint(seed,cfg)
c = cfg.model; n = c.N; m = c.M; k = c.full_neighbors;
c.ids = mod((0:n-1).' + (1:k),n) + 1;
grid = (0:m-1)'*2*pi/m; d = abs(grid-grid.'); d = min(d,2*pi-d);
W = cos(pi*(d/pi).^c.connection_exponent); c.Wfft = fft(W(:,1));
rng(seed,'twister'); x.position = rand(2,n)*c.L;
x.u = zeros(m,n); x.heading = zeros(1,n);
for i = 1:n
    x.heading(i) = 2*pi*rand; x.u(:,i) = .1*randn(m,1);
end
x.alpha = mod(grid+x.heading,2*pi);
stream = RandStream('mt19937ar','Seed',seed+100000);
order = zeros(n,k);
for i = 1:n, order(i,:) = randperm(stream,k); end
history = zeros(cfg.pre_steps+1,n,k);
history(1,:,:) = bearings(x.position,c);
for t = 1:cfg.pre_steps
    z = bearings(x.position,c);
    I = socialField(x.alpha,z,true(n,k),c,1);
    [x,~] = advance(x,I,c);
    history(t+1,:,:) = bearings(x.position,c);
end
checkpoint = x;
end

function visible = visibleMask(order,r,c)
visible = false(c.N,c.full_neighbors);
for i = 1:c.N, visible(i,order(i,1:r)) = true; end
end

function S0 = initialSocialInput(x,visible,c,kappa)
z = bearings(x.position,c);
delta = z-x.heading.';
weight = c.single_edge_weight*double(visible);
S0 = mean(sum(weight.*exp(kappa*(cos(delta)-1)),2));
end

function out = simulateBranch(checkpoint,recent,c,visible,method,params,s,diagnostics)
% Each branch updates its own 50-sample history during fully visible refresh.
x = checkpoint; n = c.N; k = c.full_neighbors; hidden = ~visible;
total = s.cycles*(s.occlusion_steps+s.refresh_steps);
series = zeros(total,10); row = 0;
for cycle = 1:s.cycles
    if any(strcmp(method,{'Hold','Linear','DampedHolt','SG','KF'}))
        forecasts = buildBearingForecasts(recent,c.dt,s.occlusion_steps,params);
        p = forecasts.(method);
    end
    for step = 1:s.occlusion_steps
        truth = bearings(x.position,c); z = truth;
        active = true(n,k); scale = 1; bearingErr = 0;
        if strcmp(method,'Drop') || strcmp(method,'StrengthOnly')
            active = visible; z(hidden) = NaN; bearingErr = NaN;
            if strcmp(method,'StrengthOnly'), scale = k/sum(visible(1,:)); end
        elseif ~strcmp(method,'Fresh')
            pred = reshape(p(step,:,:),n,k);
            z(hidden) = pred(hidden);
            error = atan2(sin(pred(hidden)-truth(hidden)),cos(pred(hidden)-truth(hidden)));
            bearingErr = mean(abs(error));
        end
        I = socialField(x.alpha,z,active,c,scale);
        if diagnostics
            oracle = socialField(x.alpha,truth,true(n,k),c,1);
            inputErr = sqrt(mean((I(:)-oracle(:)).^2));
        else
            inputErr = NaN;
        end
        [x,v] = advance(x,I,c); q = movementConsistency(x.position,v,c.L);
        row = row+1;
        series(row,:) = [q.GO,q.LO,q.C,q.movingFraction,q.meanSpeed, ...
            q.pairDistance,mean(sum(I,1)*2*pi/c.M),bearingErr,inputErr,1];
    end
    recent = zeros(s.refresh_steps,n,k);
    for step = 1:s.refresh_steps
        truth = bearings(x.position,c);
        I = socialField(x.alpha,truth,true(n,k),c,1);
        [x,v] = advance(x,I,c); q = movementConsistency(x.position,v,c.L);
        % This is the post-update observation immediately before next step.
        recent(step,:,:) = bearings(x.position,c);
        row = row+1;
        series(row,:) = [q.GO,q.LO,q.C,q.movingFraction,q.meanSpeed, ...
            q.pairDistance,mean(sum(I,1)*2*pi/c.M),NaN,0,0];
    end
    assert(max(abs(reshape(recent(end,:,:),n,k)-bearings(x.position,c)),[],'all') < 1e-12);
end
assert(row == total);
assert(sum(series(:,10)) == s.cycles*s.occlusion_steps);
out = struct('series',series,'finalState',x);
end

function z = bearings(position,c)
dx = reshape(position(1,c.ids),c.N,c.full_neighbors)-position(1,:).';
dy = reshape(position(2,c.ids),c.N,c.full_neighbors)-position(2,:).';
dx = dx-c.L*round(dx/c.L); dy = dy-c.L*round(dy/c.L);
z = mod(atan2(dy,dx),2*pi);
end

function I = socialField(alpha,z,active,c,scale)
center = z; center(~active) = 0;
delta = mod(alpha-reshape(center,1,c.N,c.full_neighbors)+pi,2*pi)-pi;
I = scale*c.single_edge_weight*sum(exp(-.5*(delta/c.kernel_sigma).^2).* ...
    reshape(active,1,c.N,c.full_neighbors),3);
end

function [x,v] = advance(x,I,c)
x.u = x.u+c.dt*(-x.u+real(ifft(fft(tanh(c.beta*x.u),[],1).*c.Wfft,[],1))/c.M- ...
    c.global_inhibition+I);
y = max(0,tanh(c.beta*x.u));
v = c.speed_scale*[sum(y.*cos(x.alpha),1);sum(y.*sin(x.alpha),1)];
moving = abs(v(1,:))+abs(v(2,:)) >= c.speed_scale*1e-9;
x.heading(moving) = mod(atan2(v(2,moving),v(1,moving)),2*pi);
x.position = mod(x.position+c.dt*v,c.L);
end

function runSmoke(seed,cfg,root)
s = cfg.schedule_candidates(1); s.cycles = 2;
[checkpoint,history,c,order] = makeCheckpoint(seed,cfg);
recent = history(end-cfg.history_steps+1:end,:,:);
full = true(c.N,c.full_neighbors);
fresh = simulateBranch(checkpoint,recent,c,full,'Fresh',struct(),s,true);
dropFull = simulateBranch(checkpoint,recent,c,full,'Drop',struct(),s,true);
assert(max(abs(fresh.series(:,1:7)-dropFull.series(:,1:7)),[],'all') < 1e-12);
visible5 = visibleMask(order,5,c); visible8 = visibleMask(order,8,c);
S5 = initialSocialInput(checkpoint,visible5,c,cfg.kappa);
S8 = initialSocialInput(checkpoint,visible8,c,cfg.kappa);
assert(S5 <= 5*c.single_edge_weight+1e-12 && S8 >= S5);
pred = simulateBranch(checkpoint,recent,c,visible5,'KF', ...
    struct('alpha',.35,'beta',.2,'phi',.85,'q',1e-5,'r',1e-4),s,true);
assert(all(isfinite(pred.series(:,1:7)),'all'));
out = fullfile(root,'results','smoke');
if ~exist(out,'dir'), mkdir(out); end
checks = struct('status','PASS','seed',seed,'steps',size(fresh.series,1), ...
    'blackoutSteps',sum(fresh.series(:,10)),'allVisibleEqualFresh',true, ...
    'S0_r5',S5,'S0_r8',S8,'ownRefreshHistory',true);
writeJson(fullfile(out,'checks.json'),checks);
end
