function track = bearingTrackStep(track,observation,dt,q,R)
% observation is empty during occlusion; wrapped innovation, Joseph covariance.
validateattributes(dt,{'numeric'},{'scalar','finite','positive'});
validateattributes(q,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(R,{'numeric'},{'scalar','finite','positive'});
if ~isempty(observation), validateattributes(observation,{'numeric'},{'scalar','finite','real'}); end
F=[1 dt;0 1]; Q=q*[dt^3/3 dt^2/2;dt^2/2 dt];
if isempty(track)
    if isempty(observation), return; end
    track.x=[mod(observation,2*pi);0]; track.P=diag([R 1]);
    track.updated=true; track.age=0; track.innovation=NaN; track.nis=NaN; return;
end
track.x=F*track.x; track.P=F*track.P*F'+Q;
if ~isfield(track,'age'), track.age=0; end
track.age=track.age+dt; track.updated=false; track.innovation=NaN; track.nis=NaN;
if ~isempty(observation)
    innovation=atan2(sin(observation-track.x(1)),cos(observation-track.x(1)));
    K=track.P(:,1)/(track.P(1,1)+R);
    track.innovation=innovation; track.nis=innovation^2/(track.P(1,1)+R);
    track.age=0; track.updated=true;
    track.x=track.x+K*innovation;
    A=eye(2)-K*[1 0]; track.P=A*track.P*A'+K*R*K';
end
track.x(1)=mod(track.x(1),2*pi);
track.P=(track.P+track.P')/2;
end
