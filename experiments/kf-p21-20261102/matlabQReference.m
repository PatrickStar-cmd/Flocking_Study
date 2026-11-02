function matlabQReference(out)
% Independent untouched MATLAB filters, actual P2 observations for every Q.
here=fileparts(mfilename('fullpath'));root=fileparts(fileparts(here));
addpath(fullfile(root,'effective-social-input'));
if nargin<1,out=fullfile(root,'effective-social-input','results','kf-p21-v1-20261102-r01','engine-tests-r01');end
path=fullfile(out,'matlab_q_reference.mat');assert(~isfile(path));
p=fullfile(root,'effective-social-input','results','kf-p2-v1-20261102-reviewed','raw','input','seed2026110101_base.mat');
a=load(p,'observed_bearing');obs=a.observed_bearing(:,[1 70 140]);
missing=false(2200,1);missing(1401:1520)=true;
refs=struct('q',{},'tau',{},'theta',{},'omega',{},'P',{},'nis',{},'age',{});
for m=1:6
 qs=[1e-6 1e-5 1e-4 1e-3 1e-2 1e-5];q=qs(m);
 if m==6,tau=10;else,tau=Inf;end
 theta=nan(2200,3);omega=theta;nis=theta;age=theta;P=nan(2200,3,2,2);
 for edge=1:3
  track=[];
  for t=1:2200
   if missing(t),y=[];else,y=obs(t,edge);end
   if isinf(tau),track=bearingTrackStep(track,y,.3,q,1e-4);
   else,track=bearingTrackStepGaussMarkov(track,y,.3,q,1e-4,tau);end
   theta(t,edge)=track.x(1);omega(t,edge)=track.x(2);P(t,edge,:,:)=track.P;
   nis(t,edge)=track.nis;age(t,edge)=track.age;
  end
 end
 refs(m)=struct('q',q,'tau',tau,'theta',theta,'omega',omega,'P',P,'nis',nis,'age',age);
end
save(path,'obs','missing','refs','-v7');fprintf('PASS reference six models,2200steps,3edges\n');
end
