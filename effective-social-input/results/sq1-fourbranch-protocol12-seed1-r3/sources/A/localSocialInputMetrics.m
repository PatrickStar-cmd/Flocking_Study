function [field,metrics] = localSocialInputMetrics(theta,weights,alpha,sigma)
%LOCALSOCIALINPUTMETRICS Evaluate one bearing-weighted social input field.
% theta and weights describe the active bearing observations for one focal
% agent. alpha is the allocentric neural grid in radians.

validateattributes(theta,{'numeric'},{'vector','real','finite'});
validateattributes(weights,{'numeric'},{'vector','real','finite','nonnegative'});
validateattributes(alpha,{'numeric'},{'vector','real','finite'});
validateattributes(sigma,{'numeric'},{'scalar','real','finite','positive'});
assert(numel(theta)==numel(weights),'theta and weights must have equal length.');

theta=mod(theta(:),2*pi);
weights=weights(:);
alpha=mod(alpha(:),2*pi);
field=zeros(size(alpha));
active=weights>0;
if any(active)
    delta=abs(alpha-theta(active).');
    delta=min(delta,2*pi-delta);
    field=exp(-0.5*(delta/sigma).^2)*weights(active);
end

mass=sum(field)*(2*pi/numel(alpha));
peak=max(field);
total=sum(weights);
if total>0
    z1=sum(weights.*exp(1i*theta))/total;
    z2=sum(weights.*exp(2i*theta))/total;
    concentration=abs(z1);
    secondConcentration=abs(z2);
    p=weights(active)/total;
    entropy=-sum(p.*log(p));
    activeAngles=sort(theta(active));
    if numel(activeAngles)<2
        maxGap=2*pi;
    else
        gaps=[diff(activeAngles); 2*pi-activeAngles(end)+activeAngles(1)];
        maxGap=max(gaps);
    end
else
    z1=0; z2=0; concentration=0; secondConcentration=0;
    entropy=0; maxGap=2*pi;
end

fieldMass=sum(field);
if fieldMass>0
    fieldZ1=sum(field.*exp(1i*alpha))/fieldMass;
    fieldZ2=sum(field.*exp(2i*alpha))/fieldMass;
else
    fieldZ1=0; fieldZ2=0;
end
metrics.totalWeight=total;
metrics.activeCount=sum(active);
metrics.inputMass=mass;
metrics.fieldPeak=peak;
metrics.z1=z1;
metrics.z2=z2;
metrics.concentration=concentration;
metrics.secondConcentration=secondConcentration;
metrics.angularEntropy=entropy;
metrics.maxAngularGap=maxGap;
metrics.angularCoverage=1-maxGap/(2*pi);
metrics.fieldConcentration=abs(fieldZ1);
metrics.fieldSecondConcentration=abs(fieldZ2);
metrics.hasDirection=total>0;
end
