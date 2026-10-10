function tr=bearingTrackBankStep(tr,obs,cfg)
% 向量化二维KF，与来源bearingTrackStep的圆周创新/Joseph更新相同。
% 接口只接受合法观测。缺失边预测，不做测量更新。
dt=cfg.dt; R=cfg.R; q=cfg.q;
assert(dt>0 && R>0 && q>=0);
old=tr.initialized & obs.selected;
tr.nis(:)=nan; tr.priorBearing(:)=nan; tr.priorVariance(:)=nan;
a=tr.P11(old)+2*dt*tr.P12(old)+dt^2*tr.P22(old)+q*dt^3/3;
b=tr.P12(old)+dt*tr.P22(old)+q*dt^2/2;
c=tr.P22(old)+q*dt;
tr.bearing(old)=mod(tr.bearing(old)+dt*tr.velocity(old),2*pi);
tr.P11(old)=a; tr.P12(old)=b; tr.P22(old)=c;
tr.priorBearing(old)=tr.bearing(old); tr.priorVariance(old)=a;
upd=old & obs.fresh;
v=atan2(sin(obs.bearing(upd)-tr.bearing(upd)),cos(obs.bearing(upd)-tr.bearing(upd)));
a=tr.P11(upd); b=tr.P12(upd); c=tr.P22(upd); s=a+R;
k1=a./s; k2=b./s;
tr.bearing(upd)=mod(tr.bearing(upd)+k1.*v,2*pi);
tr.velocity(upd)=tr.velocity(upd)+k2.*v; tr.nis(upd)=v.^2./s;
tr.P11(upd)=(1-k1).^2.*a+k1.^2*R;
tr.P12(upd)=(1-k1).*(b-k2.*a)+k1.*k2*R;
tr.P22(upd)=c-2*k2.*b+k2.^2.*(a+R);
new=~tr.initialized & obs.fresh;
tr.bearing(new)=mod(obs.bearing(new),2*pi); tr.velocity(new)=0;
tr.P11(new)=R; tr.P12(new)=0; tr.P22(new)=cfg.initialVelocityVariance;
tr.initialized(new)=true;
assert(all(isfinite(tr.bearing(tr.initialized))) && all(tr.P11(tr.initialized)>=-1e-14));
end
