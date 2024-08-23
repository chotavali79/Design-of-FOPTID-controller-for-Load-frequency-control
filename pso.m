%
% Copyright (c) 2015, Mostapha Kalami Heris & Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "LICENSE" file for license terms.
%
% Project Code: YPEA102
% Project Title: Implementation of Particle Swarm Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Cite as:
% Mostapha Kalami Heris, Particle Swarm Optimization in MATLAB (URL: https://yarpiz.com/50/ypea102-particle-swarm-optimization), Yarpiz, 2015.
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

clc;
clear;
close all;

%% Problem Definition

CostFunction = @(x) Sphere(x);        % Cost Function

nVar = 8;            % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

VarMin = [0 0 0 20 0 0 0 20];     % Lower Bound of Variables
VarMax = [5 5 5 100 5 5 5 100];

%% PSO Parameters

MaxIt = 50;      % Maximum Number of Iterations

nPop = 20;        % Population Size (Swarm Size)

% PSO Parameters
w = 1;            % Inertia Weight
wdamp = 0.99;     % Inertia Weight Damping Ratio
c1 = 2;         % Personal Learning Coefficient
c2 = 2.0;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO, 
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1 = 2.05;
% phi2 = 2.05;
% phi = phi1+phi2;
% chi = 2/(phi-2+sqrt(phi^2-4*phi));
% w = chi;          % Inertia Weight
% wdamp = 1;        % Inertia Weight Damping Ratio
% c1 = chi*phi1;    % Personal Learning Coefficient
% c2 = chi*phi2;    % Global Learning Coefficient

% Velocity Limits
for i=1:nVar
VelMax(1,i) = 0.1*(VarMax(1,i)-VarMin(1,i));
VelMin(1,i) = -VelMax(1,i);
end

%% Initialization

empty_particle.Position = [];
empty_particle.Cost = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);

GlobalBest.Cost = inf;

for i = 1:nPop
    tic
    % Initialize Position
    particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
    x= particle(i).Position;
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    %particle(i).Cost = CostFunction(particle(i).Position);
     
    Kp1=x(1,1);
    Ki1=x(1,2);
    Kd1=x(1,3);
    N1=x(1,4);
    Kp2=x(1,5);
    Ki2=x(1,6);
    Kd2=x(1,7);
    N2=x(1,8);
    
    [tout] = sim('btech_ev_1');
    particle(i).Cost = tout.e;    
    % Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest = particle(i).Best;
        
    end
    toc
    
end

BestCost = zeros(MaxIt, 1);

%% PSO Main Loop

for it = 1:MaxIt
    it
    
    for i = 1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity, VelMin);
        particle(i).Velocity = min(particle(i).Velocity, VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside = (particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside) = -particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        
        % Evaluation
        %particle(i).Cost = CostFunction(particle(i).Position);
        x= particle(i).Position;
    Kp1=x(1,1);
    Ki1=x(1,2);
    Kd1=x(1,3);
    N1=x(1,4);
    Kp2=x(1,5);
    Ki2=x(1,6);
    Kd2=x(1,7);
    N2=x(1,8);
    [tout] = sim('btech_ev_1');
    particle(i).Cost = tout.e;  
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest = particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it) = GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w = w*wdamp;
    
end

BestSol = GlobalBest;

%% Results

figure;
%plot(BestCost, 'LineWidth', 2);
semilogy(BestCost, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
