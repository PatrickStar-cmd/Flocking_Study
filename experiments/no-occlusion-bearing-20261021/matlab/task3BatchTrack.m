function [s,nis,innovation] = task3BatchTrack(s,y,fresh,dt,q,R,tau)
% Same equations and Joseph update as the scalar CV/GM implementation.
E=numel(y);
if isempty(s)
 s=struct('x1',zeros(E,1),'x2',zeros(E,1),'p11',zeros(E,1), ...
  'p12',zeros(E,1),'p22',zeros(E,1),'has',false(E,1));
end
nis=nan(E,1); innovation=nan(E,1);
if isinf(tau)
 f=dt; r=1; Q=q*[dt^3/3 dt^2/2;dt^2/2 dt];
else
 r=exp(-dt/tau); f=tau*(1-r);
 Q=q*[tau^2*dt-tau^3*(2*(1-r)-(1-r^2)/2),tau^2*(1-r)^2/2; ...
 tau^2*(1-r)^2/2,tau*(1-r^2)/2];
end
old=s.has; init=fresh & ~old;
s.x1(old)=mod(s.x1(old)+f*s.x2(old),2*pi); s.x2(old)=r*s.x2(old);
a=s.p11+2*f*s.p12+f^2*s.p22+Q(1,1);
b=r*(s.p12+f*s.p22)+Q(1,2); c=r^2*s.p22+Q(2,2);
s.p11(old)=a(old); s.p12(old)=b(old); s.p22(old)=c(old);
upd=fresh & old;
v=atan2(sin(y(upd)-s.x1(upd)),cos(y(upd)-s.x1(upd)));
den=a(upd)+R; k1=a(upd)./den; k2=b(upd)./den;
s.x1(upd)=mod(s.x1(upd)+k1.*v,2*pi); s.x2(upd)=s.x2(upd)+k2.*v;
s.p11(upd)=(1-k1).^2.*a(upd)+k1.^2*R;
s.p12(upd)=(1-k1).*(b(upd)-k2.*a(upd))+k1.*k2*R;
s.p22(upd)=c(upd)-2*k2.*b(upd)+k2.^2.*den;
innovation(upd)=v; nis(upd)=v.^2./den;
s.x1(init)=mod(y(init),2*pi); s.x2(init)=0;
s.p11(init)=R; s.p12(init)=0; s.p22(init)=1; s.has(init)=true;
end
