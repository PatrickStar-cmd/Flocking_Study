function testBearingTrackGaussMarkov()
dt=.3; q=1e-4; R=1e-4;

% Infinite correlation time must recover the existing constant-velocity KF.
a=[]; b=[];
observations={6.27,.02,.04,[],[],.08};
for i=1:numel(observations)
    a=bearingTrackStep(a,observations{i},dt,q,R);
    b=bearingTrackStepGaussMarkov(b,observations{i},dt,q,R,Inf);
    assert(isequaln(a,b));
end

% Missing observations must apply the exact Gauss-Markov prior only.
tau=10;
t=bearingTrackStepGaussMarkov([],6.27,dt,q,R,tau);
t=bearingTrackStepGaussMarkov(t,.02,dt,q,R,tau);
rho=exp(-dt/tau);
F=[1 tau*(1-rho);0 rho];
Q=q*[tau^2*dt-tau^3*(2*(1-rho)-(1-rho^2)/2), ...
    tau^2*(1-rho)^2/2; ...
    tau^2*(1-rho)^2/2, tau*(1-rho^2)/2];
expectedP=F*t.P*F'+Q;
expectedX=F*t.x;
p=bearingTrackStepGaussMarkov(t,[],dt,q,R,tau);
assert(max(abs(p.P-expectedP),[],'all')<1e-12);
assert(abs(atan2(sin(p.x(1)-expectedX(1)),cos(p.x(1)-expectedX(1))))<1e-12);
assert(abs(p.x(2)-expectedX(2))<1e-12);
assert(~p.updated && isnan(p.nis) && p.age==dt);
assert(all(eig(p.P)>-1e-12));

u=bearingTrackStepGaussMarkov(p,.03,dt,q,R,tau);
assert(u.updated && u.age==0 && isfinite(u.nis));
assert(all(eig(u.P)>-1e-12));
disp('PASS: Gauss-Markov bearing model contract');
end
