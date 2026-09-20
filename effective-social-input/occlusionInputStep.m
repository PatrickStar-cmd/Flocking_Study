function [input,state,diag] = occlusionInputStep(observation,alpha,state,opt)
% Only sanitized observations enter here: hidden bearing/weight MUST be NaN.
n=size(observation.fresh,1); ns=size(alpha,2);
if isempty(state)
 state.tracks=cell(n); state.bearing=nan(n); state.weight=nan(n); state.age=inf(n);
end
input=zeros(ns,n); diag.active=zeros(n); diag.variance=nan(n); diag.center=nan(n);
for i=1:n
 for j=find(observation.selected(i,:))
  fresh=observation.fresh(i,j);
  if fresh
   y=observation.bearing(i,j); w=observation.weight(i,j);
   assert(isfinite(y)&&isfinite(w)); state.bearing(i,j)=y; state.weight(i,j)=w; state.age(i,j)=0;
  else
   assert(isnan(observation.bearing(i,j))&&isnan(observation.weight(i,j)));
   y=[]; state.age(i,j)=state.age(i,j)+opt.dt;
  end
  if strcmp(opt.mode,'kalman-gm')
   track=bearingTrackStepGaussMarkov(state.tracks{i,j},y,opt.dt,opt.gmQ,opt.R,opt.velocityTau);
  else
   track=bearingTrackStep(state.tracks{i,j},y,opt.dt,opt.q,opt.R);
  end
  state.tracks{i,j}=track;
  if isempty(track), continue; end
  if fresh
   center=y; variance=0;
  elseif strcmp(opt.mode,'none') || state.age(i,j)>opt.maxAge
   continue;
  elseif strcmp(opt.mode,'hold')
   center=state.bearing(i,j); variance=0;
  elseif strcmp(opt.mode,'mean') || strcmp(opt.mode,'kf-mean')
   if track.P(1,1)>opt.maxVariance, continue; end
   center=track.x(1); variance=0;
  elseif strcmp(opt.mode,'kalman') || strcmp(opt.mode,'kalman-gm')
   if track.P(1,1)>opt.maxVariance, continue; end
   center=track.x(1); variance=track.P(1,1);
  else
   error('Unknown prediction mode: %s',opt.mode);
  end
  diag.center(i,j)=center;
  delta=atan2(sin(alpha(i,:)-center),cos(alpha(i,:)-center));
  kernel=exp(-.5*(delta/opt.sigma).^2);
  if variance>0
   % Circular heat convolution: preserves discrete mass, tends to uniform.
   frequency=[0:floor(ns/2),-ceil(ns/2)+1:-1];
   kernel=real(ifft(fft(kernel).*exp(-.5*variance*frequency.^2)));
   kernel=max(0,kernel);
  end
  input(:,i)=input(:,i)+opt.h0*state.weight(i,j)*kernel';
  diag.active(i,j)=1; diag.variance(i,j)=track.P(1,1);
 end
end
diag.age=state.age; diag.weight=state.weight;
end
