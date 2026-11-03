function m=movementConsistency(position,velocity,L)
% Paper GO/LO, with self plus five nearest neighbors for LO.
n=size(position,2); assert(isequal(size(position),[2 n]) && isequal(size(velocity),[2 n]));
speed=hypot(velocity(1,:),velocity(2,:)); unit=zeros(2,n); moving=speed>1e-12;
unit(:,moving)=velocity(:,moving)./speed(moving);
dx=position(1,:).'-position(1,:); dy=position(2,:).'-position(2,:);
dx=dx-L*round(dx/L); dy=dy-L*round(dy/L); distance=hypot(dx,dy);
m.GO=norm(sum(unit,2))/n;
otherDistance=distance; otherDistance(1:n+1:end)=inf;
% Stable column sorting breaks exact ties by agent ID. Add self explicitly.
[~,order]=sort(otherDistance,2,'ascend'); k=min(5,n-1); ids=order(:,1:k);
cx=unit(1,:).'+sum(reshape(unit(1,ids),n,k),2);
cy=unit(2,:).'+sum(reshape(unit(2,ids),n,k),2);
m.LO=mean(hypot(cx,cy))/(k+1); m.C=min(m.GO,m.LO); m.movingFraction=mean(moving); m.meanSpeed=mean(speed);
m.pairDistance=mean(distance(triu(true(n),1)))/L;
end
