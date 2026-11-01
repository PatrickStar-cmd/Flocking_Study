function resultsDir=runKFDiagnosticBatch(seeds,suffix)
% Reviewed runner: every condition saved and audited; no failed batch success.
if nargin<1, seeds=2026110101:2026110105; end
if nargin<2, suffix=char(datetime('now','Format','yyyyMMdd-HHmmss')); end
resultsDir=runReviewedBatch(seeds,suffix);
end
