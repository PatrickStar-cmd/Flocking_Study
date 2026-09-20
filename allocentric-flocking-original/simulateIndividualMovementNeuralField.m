
%SIMULATEINDIVIDUALMOVEMENT  Single‐agent ring‐attractor (optional target)
%   If targetInitialPos is empty runs free motion.

  
% Author: Mohammad Salahshour


function [xPos, yPos, TX, TY, uArr, headings] = simulateIndividualMovementNeuralField( ...
    N_s, T, dt, v0, L, W, h_b, beta, h0, sigma, ...
    initialPos, initialHeading, targetSpeed, targetInitialPos, periodicFlag, allocentricFlag)

  hasTarget = ~isempty(targetInitialPos) ;

  %% Preallocate
  xPos     = zeros(1,T+1);
  yPos     = zeros(1,T+1);
  headings = zeros(1,T+1);
  uArr     = zeros(N_s,T+1);

  if hasTarget
    TX = zeros(1,T+1);
    TY = zeros(1,T+1);
  else
    TX = []; TY = [];
  end

  % initialize
  xPos(1)     = initialPos(1);
  yPos(1)     = initialPos(2);
  headings(1) = initialHeading;
  uArr(:,1)   = 0.1*randn(N_s,1);
  if hasTarget
    TX(1) = targetInitialPos(1);
    TY(1) = targetInitialPos(2);
  end

  % ring angles (column)
  alpha0 = linspace(0,2*pi,N_s+1); alpha0=alpha0(1:end-1).';
  if allocentricFlag
    alpha = alpha0;
  else
    alpha = mod(alpha0 + initialHeading,2*pi);
  end

  f = @(u)tanh(beta*u);  % activation

  %% Time‐loop
  for t = 1:T
    % external input
    if hasTarget
      dx = TX(t) - xPos(t);
      dy = TY(t) - yPos(t);
      if periodicFlag
        if abs(dx)>L/2, dx = dx - sign(dx)*L; end
        if abs(dy)>L/2, dy = dy - sign(dy)*L; end
      end
      theta = atan2(dy,dx);
      if theta<0, theta=theta+2*pi; end
      dAngs = abs(alpha - theta);
      dAngs = min(dAngs,2*pi-dAngs);
      Iext  = h0 * exp(-0.5*(dAngs.^2)/(sigma^2));
    else
      Iext = zeros(N_s,1);
    end

    % Euler update ring
    u    = uArr(:,t);
    net  = (W*f(u))/N_s;
    du   = -u + net - h_b + Iext;
    uNew = u + dt*du;
    uArr(:,t+1) = uNew;

    % read‐out heading
    aNew = f(uNew); aNew(aNew<0)=0;
    cx = sum(aNew .* cos(alpha));
    cy = sum(aNew .* sin(alpha));
    if abs(cx)<1e-9 && abs(cy)<1e-9
      newH = headings(t);
    else
      newH = atan2(cy,cx);
      if newH<0, newH=newH+2*pi; end
    end
    headings(t+1)=newH;

    % rotate ring if egocentric
    if ~allocentricFlag
      alpha = mod(alpha - headings(t) + newH, 2*pi);
    end

    % move agent
    xPos(t+1)=xPos(t)+dt*v0*cx;
    yPos(t+1)=yPos(t)+dt*v0*cy;
    if periodicFlag
      xPos(t+1)=mod(xPos(t+1),L);
      yPos(t+1)=mod(yPos(t+1),L);
    else
      xPos(t+1)=min(max(xPos(t+1),0),L);
      yPos(t+1)=min(max(yPos(t+1),0),L);
    end

    % move target if any
    if hasTarget
      if rand<0.5
        TX(t+1)=mod(TX(t)+targetSpeed*sign(randn),L);
        TY(t+1)=TY(t);
      else
        TY(t+1)=mod(TY(t)+targetSpeed*sign(randn),L);
        TX(t+1)=TX(t);
      end
    end
        %% Visualization 
    if mod(t,10)==0
        figure(1)
        figure(1); clf; hold on;
        % Plot agents
       
            scatter(xPos(:,t+1), yPos(:,t+1), 80, 'filled','MarkerFaceColor',[0,0.2,0.8]);
        
        % Plot targets
        if ~isempty(TX)
            scatter(TX(:,t+1), TY(:,t+1), 80, 's','filled','MarkerFaceColor',[0.8,0,0.2]);
        end
        axis([0 L 0 L]); axis square;
        title(sprintf('Time step %d (beta=%.2f)', t, beta));
        drawnow;



    end

  end
end
