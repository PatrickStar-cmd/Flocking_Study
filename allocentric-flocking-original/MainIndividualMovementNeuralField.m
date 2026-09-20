% MAININDIVIDUALMOVEMENT 

clear; clc;

%% 1) Parameters
N_s          = 100;           % ring neurons
T            = 1000;          % time steps
dt           = 0.3;             % step size
v0           = 0.05;            % agent speed constant
L            = 100;          % arena size
h_b          = 0;             % inhibition
beta         = 20;           % inverse‐noise
h0           = 0.25;        % target field strength
sigma        = 0.4;      % receptive field width
periodicFlag = true;          % whether the space is periodic or not.
allocFlag    = true;          % allocentric if true and egocentric is false

% connectivity (cosine kernel)
nu = 0.5;
a_tmp = linspace(0,2*pi,N_s+1);
a_tmp = a_tmp(1:end-1);
W = zeros(N_s);
for i=1:N_s
  for j=1:N_s
    d = abs(a_tmp(i)-a_tmp(j)); d = min(d,2*pi-d);
    W(i,j) = cos((d/pi)^nu*pi);
  end
end

% initial agent
initPos     = [500;500];
initHead    = pi;  

% choose target
targetSpeed = 0;      
targetInit  = [50;50];

%% 2) Run simulation
[xPos,yPos,TX,TY,uArr,head] = simulateIndividualMovementNeuralField( ...
  N_s,T,dt,v0,L,W,h_b,beta,h0,sigma,...
  initPos,initHead,targetSpeed,targetInit,...
  periodicFlag,allocFlag);

hasTarget = ~isempty(TX);

save('SamplerunIndividualMovement')