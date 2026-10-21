function [e,isOutlier,z,state]=drawNonGaussianBearingNoise(state,cfg,shape)
% Simulator-only generator: unit population variance, zero mean, heavy tails.
% Gaussian base draws match the prior pilot; mixture labels use another stream.
p=cfg.nonGaussianOutlierProbability; ratio=cfg.nonGaussianOutlierScale;
assert(p>=0 && p<1 && ratio>=1);
if ~isfield(state,'noiseStream')
    state.noiseStream=RandStream('mt19937ar','Seed',cfg.seed+15000);
    state.initialNoiseStreamState=state.noiseStream.State;
    state.mixtureStream=RandStream('mt19937ar','Seed',cfg.seed+25000);
    state.initialMixtureStreamState=state.mixtureStream.State;
end
z=randn(state.noiseStream,shape);
isOutlier=rand(state.mixtureStream,shape)<p;
scale=ones(shape); scale(isOutlier)=ratio;
e=z.*scale/sqrt((1-p)+p*ratio^2);
end
