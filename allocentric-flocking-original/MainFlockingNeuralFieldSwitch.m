


% MAINCOLLECTIVEFLOCKING  Example: allocentric flocking with random allo-ego switch



% base params
N_agents       = 20;
N_s            = 100;
T              = 8000;
dt             = 0.3;
v0             = 0.05;
L              = 1000;
h_b            = 0;
beta           = 1000;
h_t            = 0.24;          % total social attraction
h0s            = h_t / N_agents;
sigma          = 0.4;
egocentricProb= 0.5;            

% connectivity W
nu = 0.5;
alpha0 = linspace(0,2*pi,N_s+1); alpha0=alpha0(1:end-1);
W = zeros(N_s,N_s);
for i=1:N_s
  for j=1:N_s
    d = abs(alpha0(i)-alpha0(j));
    d = min(d,2*pi-d);
    W(i,j) = cos((d/pi)^nu * pi);
  end
end

% random initial positions
initialXY = rand(2,N_agents)*L;

% run
[xPos,yPos] = simulateCollectiveFlockingNeuralFieldAlloEgoSwitch( ...
  N_agents, N_s, T, dt, v0, L, W, h_b, beta, ...
  h0s, sigma, initialXY, egocentricProb);


save('SamplerunCollectiveMovementSwitch')