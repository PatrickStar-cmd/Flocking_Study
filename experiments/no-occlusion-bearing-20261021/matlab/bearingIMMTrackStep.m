function [s,d] = bearingIMMTrackStep(s,y,fresh,dt,qPair,R,transition)
% Two angular-CV models with different process-noise intensities.
% Only observations and freshness enter this function. No truth/range input.
% transition(i,j) is a probability per MODEL STEP, not per second.
% Circular local-Gaussian approximation: shortest residuals and tangent
% covariances, not an exact Bayesian filter for an arbitrary circular density.
E=numel(y); y=y(:); fresh=fresh(:);
assert(numel(qPair)==2 && all(qPair>0) && R>0 && dt>0);
assert(isequal(size(transition),[2,2]) && all(transition(:)>=0));
assert(max(abs(sum(transition,2)-1))<1e-12);
assert(all(isfinite(y(fresh))) && all(isnan(y(~fresh))));
if isempty(s)
 z=struct('x1',zeros(E,1),'x2',zeros(E,1),'p11',zeros(E,1), ...
  'p12',zeros(E,1),'p22',zeros(E,1),'has',false(E,1));
 s.models={z,z}; s.mu=repmat([.5,.5],E,1);
end
old=s.models{1}.has; init=fresh & ~old;
priorMu=s.mu*transition; mixed=s.models; weakMix=0;
for j=1:2
 w=s.mu.*transition(:,j).'; w=w./priorMu(:,j);
 [theta,weak]=circleMean(s.models{1}.x1,s.models{2}.x1,w);
 omega=w(:,1).*s.models{1}.x2+w(:,2).*s.models{2}.x2;
 p11=zeros(E,1); p12=p11; p22=p11;
 for i=1:2
  a=atan2(sin(s.models{i}.x1-theta),cos(s.models{i}.x1-theta));
  b=s.models{i}.x2-omega;
  p11=p11+w(:,i).*(s.models{i}.p11+a.^2);
  p12=p12+w(:,i).*(s.models{i}.p12+a.*b);
  p22=p22+w(:,i).*(s.models{i}.p22+b.^2);
 end
 mixed{j}.x1(old)=theta(old); mixed{j}.x2(old)=omega(old);
 mixed{j}.p11(old)=p11(old); mixed{j}.p12(old)=p12(old);
 mixed{j}.p22(old)=p22(old); weakMix=weakMix+sum(weak & old);
end
innov=nan(E,2); variance=nan(E,2); nis=nan(E,2);
for j=1:2
 variance(:,j)=mixed{j}.p11+2*dt*mixed{j}.p12+ ...
  dt^2*mixed{j}.p22+qPair(j)*dt^3/3+R;
 [s.models{j},nis(:,j),innov(:,j)]=task3BatchTrack( ...
  mixed{j},y,fresh,dt,qPair(j),R,Inf);
end
s.mu(old,:)=priorMu(old,:);
upd=fresh & old;
% Stable log likelihood prevents underflow during a near-pi bearing jump.
logWeight=log(priorMu(upd,:))-.5*(log(2*pi*variance(upd,:))+ ...
 innov(upd,:).^2./variance(upd,:));
if any(upd)
 weight=exp(logWeight-max(logWeight,[],2));
 s.mu(upd,:)=weight./sum(weight,2);
end
s.mu(init,:)=repmat([.5,.5],sum(init),1);
[theta,weak]=circleMean(s.models{1}.x1,s.models{2}.x1,s.mu);
omega=s.mu(:,1).*s.models{1}.x2+s.mu(:,2).*s.models{2}.x2;
p11=zeros(E,1); p12=p11; p22=p11;
for i=1:2
 a=atan2(sin(s.models{i}.x1-theta),cos(s.models{i}.x1-theta));
 b=s.models{i}.x2-omega;
 p11=p11+s.mu(:,i).*(s.models{i}.p11+a.^2);
 p12=p12+s.mu(:,i).*(s.models{i}.p12+a.*b);
 p22=p22+s.mu(:,i).*(s.models{i}.p22+b.^2);
end
s.x1=theta; s.x2=omega; s.p11=p11; s.p12=p12; s.p22=p22;
s.has=s.models{1}.has;
d=struct('innovation',innov,'innovationVariance',variance,'nis',nis, ...
 'weakCircularMeans',weakMix+sum(weak & s.has));
end

function [theta,weak]=circleMean(a,b,w)
cx=w(:,1).*cos(a)+w(:,2).*cos(b);
cy=w(:,1).*sin(a)+w(:,2).*sin(b);
weak=hypot(cx,cy)<1e-10;
theta=mod(atan2(cy,cx),2*pi);
% An almost-antipodal mixture has no well-defined circular mean. Retain the
% more probable component, with model 1 winning an exact tie; count it.
theta(weak)=a(weak); pick=weak & w(:,2)>w(:,1); theta(pick)=b(pick);
end
