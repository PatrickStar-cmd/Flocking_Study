function matlabParityReference(out)
% Independent references call unchanged stock MATLAB filters, not the port.
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
addpath(fullfile(root,'effective-social-input'));
addpath(fullfile(root,'experiments','kf-diagnostic-20261101'));
if nargin<1, out=fullfile(root,'effective-social-input','results','kf-p2-v1-20261102-reviewed'); end
path=fullfile(out,'raw','matlab_reference.mat'); assert(~isfile(path));
d=jsondecode(fileread(fullfile(out,'design.json')));
refs=struct('case_id',{},'D',{},'edge_indices',{},'cv',{},'gm',{});
first=min(d.seeds); edge_indices=[1 70 140];
for c=1:numel(d.cases)
 ca=d.cases(c); if ca.seed~=first, continue; end
 a=load(fullfile(out,ca.source),'observed_bearing'); obs=a.observed_bearing(:,edge_indices);
 for D=d.durations(:).'
  masks=false(size(obs,1),1); masks(ca.t0:ca.t0+D-1)=true;
  ref=struct('case_id',ca.case_id,'D',D,'edge_indices',edge_indices, ...
   'cv',stock(obs,masks,d,false),'gm',stock(obs,masks,d,true));
  refs(end+1)=ref; %#ok<AGROW>
 end
end
save(path,'refs','-v7'); fprintf('PASS stock MATLAB references %d cases x2 filters x3 edges\n',numel(refs));
end

function data=stock(obs,missing,d,gm)
[T,E]=size(obs); keys={'theta','p11','p12','p22','nis','age'}; data=struct();
for k=1:numel(keys), data.(keys{k})=nan(T,E); end
for e=1:E
 track=[];
 for t=1:T
  if missing(t), y=[]; else, y=obs(t,e); end
  if gm
   track=bearingTrackStepGaussMarkov(track,y,d.dt,d.q_gm,d.R,d.tau);
  else
   track=bearingTrackStep(track,y,d.dt,d.q_cv,d.R);
  end
  data.theta(t,e)=track.x(1); data.p11(t,e)=track.P(1,1);
  data.p12(t,e)=track.P(1,2); data.p22(t,e)=track.P(2,2);
  data.nis(t,e)=track.nis; data.age(t,e)=track.age;
 end
end
end
