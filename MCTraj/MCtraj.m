% clear all
clearvars
clearvars -GLOBAL
close all
% set(0,'DefaultFigureWindowStyle','docked')
global C

addpath ../geom2d/geom2d

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665                       %metres (32.1740 ft) per s²

nTraj = 10;
nSims = 100;

InitalAngle = 55 * pi / 180;
PlotTraj = 1;
MaxC = 100;
doPlot = 1;
V0 = 1000;
g = 1;
dt = 1; 

Wind = @UniformRandWind;
WindParas = [10];

xl = [];
  
    for n = 1: nSims
        x(1, :) = zeros(1, nTraj);
        y(1, :) = zeros(1, nTraj);
        time(1,:)  = zeros(1, nTraj);
        Vx(1:nTraj) = V0 * cos(InitalAngle);
        Vy(1:nTraj) = V0 * sin(InitalAngle);

        for c=2:MaxC

            dvx = Wind(nTraj,WindParas)*dt;
            Vx = Vx + dvx;
            dx = Vx * dt;
            dy = Vy * dt; 

        x(c,:) = x(c - 1,:) + dx;
        y(c,:) = y(c - 1,:) + dy;
        time(c,:) = dt + time(c-1,:); 
        
        if max(y(c)) < 0
            break
            dvy = -g * dt;
            Vy = Vy + dvy;
            dy = Vy * dt + g * dt^2 / 2;
            
            if max(y(c - 1,:) + dy) <= 0
                break
            end
            
            x(c,:) = x(c - 1,:) + dx;
            y(c,:) = y(c - 1,:) + dy;

        end

        xl = [xl x(end, :)];
    end
end
    
%plot of the trajectories position over time 
magnitudePOS = sqrt(x.^2 + y.^2);
figure (1)
plot(time, magnitudePOS)
title('Plot of the Position of the Trajectories Over Time')
xlabel ('Time (seconds)')
ylabel ('Position (meters)')

%plot of the trajectories speed 
%through each iteration above, there is the addition of a randomly
%distrubited change in velccity which is added to the current wind
%velcocity. The changes of these velcities through each step can be seen
%below
magnitudeVEL = sqrt(Vx.^2 + Vy.^2);
figure (2)
plot(magnitudeVEL)
title('Plot of the Velocities Randomly Distributed Over Ten Iterations')
xlabel ('Ten Different Steps of Changing Wind Patterns')
ylabel ('Wind Speed (m/s)')

