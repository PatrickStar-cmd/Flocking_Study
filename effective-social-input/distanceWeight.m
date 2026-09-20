function [weights,normalized] = distanceWeight(distances,cfg)
%DISTANCEWEIGHT Candidate bounded exponential range weight.
% Distances must be observable estimates, never hidden true states in occlusion.
validateattributes(distances,{'numeric'},{'real','finite','nonnegative'});
if nargin<2 || ~isfield(cfg,'distanceWeightMode'), cfg.distanceWeightMode='equal'; end
if strcmp(cfg.distanceWeightMode,'equal')
    weights=ones(size(distances));
elseif strcmp(cfg.distanceWeightMode,'exponential')
    if ~isfield(cfg,'distanceWeightScaleFraction') || cfg.distanceWeightScaleFraction<=0
        error('distanceWeightScaleFraction must be positive.');
    end
    scale=cfg.distanceWeightScaleFraction*cfg.arenaSize;
    weights=exp(-distances/scale);
else
    error('Unknown distanceWeightMode.');
end
if isempty(distances)
    normalized=weights;
elseif strcmp(cfg.distanceWeightMode,'exponential')
    stable=exp(-(distances-min(distances))/scale);
    normalized=stable/sum(stable);
else
    normalized=weights/sum(weights);
end
end
