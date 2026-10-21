function audit=testNoOcclusionHandoff(outDir)
% Run from an isolated MATLAB path. This is software evidence only.
assert(~isfolder(outDir));
here=fileparts(mfilename('fullpath')); addpath(here,'-begin');
audit=struct();
audit.kf=testNoisyBearingKFPilot();
audit.nonGaussian=testNonGaussianBearingKFPilot();
audit.immTrack=testBearingIMMTrack();
audit.closedLoop=testBearingIMMClosedLoop();
rows=runNoOcclusionHandoff(outDir,'smoke');
assert(numel(rows)==16);
c=jsondecode(fileread(fullfile(fileparts(here),'configs','base-config.json')));
c.seed=6299; c.NSteps=80; c.burnIn=40; c.pilotSampleStride=20;
c.bearingMeasurementStd=0; c.noisyBearingMode='raw';
a=simulateOracleBearingDiagnostic(c); b=simulateBearingIMMClosedLoopPilot(c);
assert(isequal(a.trajectory,b.trajectory));
audit.oracleEqualsZeroNoiseRaw=true;
audit.status='PASS'; audit.evidence='software checks only';
fid=fopen(fullfile(outDir,'software-test-audit.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(audit,PrettyPrint=true));
end
