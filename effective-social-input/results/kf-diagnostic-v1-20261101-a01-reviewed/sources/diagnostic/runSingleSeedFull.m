function resultsDir=runSingleSeedFull(seed,suffix)
% Nine paired conditions for one developmental seed, using the reviewed runner.
if nargin<1, seed=2026110101; end
if nargin<2, suffix=sprintf('single-%d-%s',seed,char(datetime('now','Format','yyyyMMdd-HHmmss'))); end
resultsDir=runReviewedBatch(seed,suffix);
end
