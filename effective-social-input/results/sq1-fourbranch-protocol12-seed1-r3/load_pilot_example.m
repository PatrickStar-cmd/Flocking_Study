%LOAD_PILOT_EXAMPLE Locate a completed r3 run without absolute path assumptions.
here=fileparts(mfilename('fullpath'));
if isfile(fullfile(here,'fresh_seed1.mat'))
    resultDir=here;
else
    listing=dir(fullfile(here,'results','*','metadata.json')); resultDir='';
    for idx=numel(listing):-1:1
        candidate=listing(idx).folder;
        test=jsondecode(fileread(fullfile(candidate,'tests.json')));
        if all(strcmp({test.status},'PASS')), resultDir=candidate; break; end
    end
    assert(~isempty(resultDir),'No completed all-PASS result directory found.');
end
x=load(fullfile(resultDir,'hold_seed1.mat'),'run'); run=x.run;
assert(isscalar(run.checkpoint6000)&&isscalar(run.checkpoint6200));
disp(readtable(fullfile(resultDir,'summary.csv')));
fprintf('Loaded hold: %d steps, %d state frames; cp diagnostics %d agents\n',...
    numel(run.segment.data.step),size(run.segment.state.position,3),numel(run.checkpoint6200.maskDiagnostics));
% If the source package is on the MATLAB path, this checks public resume:
if exist('resumeProtocolFromCheckpoint','file')==2
    resumed=resumeProtocolFromCheckpoint(run.checkpoint6200,run.config,6600);
    difference=max(abs(resumed.finalNeuralState-run.segment.finalNeuralState),[],'all');
    fprintf('cp6200 -> 6600 final neural-state difference: %.17g\n',difference);
    assert(difference==0);
end
