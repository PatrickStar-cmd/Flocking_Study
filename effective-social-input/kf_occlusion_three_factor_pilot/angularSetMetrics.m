function s=angularSetMetrics(theta)
theta=theta(:); s.count=numel(theta); s.R1=0; s.R2=0;
s.direction=nan; s.direction_defined=false; s.maxGap=nan;
if isempty(theta), return; end
z=mean(exp(1i*theta)); s.R1=abs(z); s.R2=abs(mean(exp(2i*theta)));
s.direction_defined=s.R1>1e-8;
if s.direction_defined, s.direction=angle(z); end
a=sort(mod(theta,2*pi)); s.maxGap=max([diff(a);2*pi-a(end)+a(1)]);
end
