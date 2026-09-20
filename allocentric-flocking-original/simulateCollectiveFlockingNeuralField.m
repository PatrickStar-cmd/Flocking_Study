
% SIMULATECOLLECTIVEFLOCKING  ring‐attractor flocking (allocentric/egocentric)
%

% Author: Mohammad Salahshour

% [xPos,yPos,uArr,headings] = simulateCollectiveFlocking( ...
%     N_agents,N_s,T,dt,v0,L,W,h_b,beta,h0s,sigma, ...
%     initialXY,periodicFlag,allocentricFlag)
%
% Inputs:
%   N_agents       Number of agents
%   N_s            Number of neurons per agent
%   T, dt, v0, L   Simulation time, integration constant, speed constant,
%   and arena size.
%   W  (N_s×N_s)   ring connectivity
%   h_b, beta      inhibition & inverse noise
%   h0s            social‐attraction
%   sigma          receptive‐field width
%   initialXY(2×N) initial agent positions
%   periodicFlag   wrapping on/off
%   allocentricFlag whether the agents have an allocentric (1) or
%   egocentric (0) perception of space.
% 
 function [xPos,yPos,uArr,headings] = simulateCollectiveFlockingNeuralField( ...
    N_agents, N_s, T, dt, v0, L, W, h_b, beta, h0s, sigma, ...
    initialXY, allocentricFlag)


%% allocate
xPos     = zeros(N_agents,T+1);
yPos     = zeros(N_agents,T+1);
headings = zeros(N_agents,T+1);
uArr     = zeros(N_s,N_agents,T+1);

% init
for a=1:N_agents
  xPos(a,1)=initialXY(1,a);
  yPos(a,1)=initialXY(2,a);
  headings(a,1)=2*pi*rand;
  uArr(:,a,1)=0.1*randn(N_s,1);
end

% ring angles
alpha0 = linspace(0,2*pi,N_s+1); alpha0=alpha0(1:end-1);
alpha = repmat(alpha0,N_agents,1);

  alpha = mod(alpha + headings(:,1),2*pi);


f = @(u)tanh(beta*u);

%% time‐loop
for t=1:T
  Iext = zeros(N_s,N_agents);

  % social input
  for a=1:N_agents
    xa=xPos(a,t); ya=yPos(a,t);
    for b=1:N_agents
      if b==a, continue; end
      xb=xPos(b,t); yb=yPos(b,t);
      dx=xb-xa; dy=yb-ya;
      if abs(dx)>L/2, dx=dx-sign(dx)*L; end
      if abs(dy)>L/2, dy=dy-sign(dy)*L; end
      theta=atan2(dy,dx); if theta<0, theta=theta+2*pi; end
      dAngs = abs(alpha(a,:)-theta);
      dAngs = min(dAngs,2*pi-dAngs);
      Iext(:,a) = Iext(:,a) + ...
        h0s*exp(-0.5*(dAngs.^2)/sigma^2).';
    end
  end

  % update rings
  for a=1:N_agents
    u = uArr(:,a,t);
    aout = f(u);
    uArr(:,a,t+1) = u + dt*( -u + (W*aout)/N_s - h_b + Iext(:,a) );
  end

  % headings & optional egocentric‐ring rotate
  for a=1:N_agents
    u = uArr(:,a,t+1);
    aout = f(u); aout(aout<0)=0;
    cx = sum(aout.*cos(alpha(a,:).'));
    cy = sum(aout.*sin(alpha(a,:).'));
    if abs(cx)+abs(cy)<1e-9
      newH = headings(a,t);
    else
      newH = atan2(cy,cx);
      if newH<0, newH=newH+2*pi; end
    end
    headings(a,t+1)=newH;
    if ~allocentricFlag
      alpha(a,:) = mod(alpha(a,:) ...
                    -headings(a,t)+newH,2*pi);
    end
  end

  % move agents
  for a=1:N_agents
    u   = uArr(:,a,t+1);
    aout= f(u); aout(aout<0)=0;
    cx  = sum(aout.*cos(alpha(a,:).'));
    cy  = sum(aout.*sin(alpha(a,:).'));
    xPos(a,t+1)=xPos(a,t)+dt*v0*cx;
    yPos(a,t+1)=yPos(a,t)+dt*v0*cy;
    xPos(a,t+1)=mod(xPos(a,t+1),L);
    yPos(a,t+1)=mod(yPos(a,t+1),L);

  end
      %% Visualization 
    if mod(t,10)==0
        figure(1)
        figure(1); clf; hold on;
        % Plot agents
        if N_agents>0
            scatter(xPos(:,t+1), yPos(:,t+1), 80, 'filled','MarkerFaceColor',[0,0.2,0.8]);
        end

        axis([0 L 0 L]); axis square;
       %axis([min(xPos(:,t+1)) max(xPos(:,t+1)) min(yPos(:,t+1)) max(yPos(:,t+1))]); axis square;
        title(sprintf('Time step %d (beta=%.2f)', t, beta));
        drawnow;



    end

end
end

