%Thomas Kantner 104577846
%Rocket Trajectory Simulation

clear; close all; clc;

%Physical Parameters
g0 = 9.81; % [m/s^2]

%Inputs
m_dry = (322 - 65)/1000; % Dry mass of rocket [g]
rail_height = 4; % Launch Rail height [ft]
theta = 0; % Launch Angle [deg]
theta = theta*pi/180; % Launch Angle [rad]
d_body = 5.4/100; % Body Diameter [m];
A_body = d_body^2*pi/4; % Body Area [m^2]
d_parachute = 12; % Parachute Diameter [in]
d_parachute = d_parachute*.0254; % Parachute Diameter [m]
A_parachute = d_parachute^2*pi/4; % Parachute Area [m^2]
CD = 0.7; % Rocket Drag coeff [-]
CD_para = 1.5; % Parachute Drag coeff [-]

%Motor Specs
%Source: http://www.nar.org/SandT/pdf/Aerotech/F32T.pdf
I_tot = 56.9; % Total impulse [N*s]
m_prop = 25.8/1000; % Propellent mass [kg]
m_motor = (65 - m_prop)/1000; % Mass of motor after firing [kg]
delay = 8; % Ejection Delay [s]
T_avg = 34.1; % Avg thrust [N]
THRUST = [
0.0 0.00669778
0.17 0.236738
0.192 0.794052
0.202 2.68628
0.204 3.57052
0.208 7.96887
0.21 11.79011
0.212 16.82641
0.22 39.2314
0.224 46.779
0.228 51.8968
0.232 55.296
0.236 56.9612
0.258 58.4048
0.288 55.2426
0.316 53.4728
0.394 49.5892
0.458 47.4002
0.73 41.2237
0.904 39.3815
1.11 35.6689
1.198 35.0946
1.22 33.8331
1.2599 33.5013
1.4779 27.3792
1.4899 27.7527
1.5099 25.4378
1.6759 11.26592
1.7679 5.94688
1.8519 2.91383
1.9519 0.745781
2.0639 0.0  ];
    
%Intermediate Calculations
w_prop = m_prop*g0; % Propellent weight [N]
Isp_avg = I_tot/w_prop; % [s]
m_dot_avg = T_avg/(Isp_avg*g0); % Avg mass flow rate [kg/s]

dt = 0.001; % [s]
totTime = 80; % [s]
numSteps = ceil(totTime/dt);

%Allocate Vectors
TIME = zeros(1, numSteps); % [s]
ACCELERATION = zeros(1, numSteps); % [m/s^2]
VELOCITY = zeros(1, numSteps); % [m/s]
ALTITUDE = zeros(1, numSteps); % [m]
MASS = zeros(1, numSteps); % [kg]
F_NET = zeros(1, numSteps); % [N]

%Initial Conditions
t0 = 0; % [s]
a0 = 0; % [m/s^2]
v0 = 0; % [m/s]
y0 = 0; % [m]
T0 = 0; % [N]
TIME(1) = t0;
ACCELERATION(1) = a0;
VELOCITY(1) = v0;
ALTITUDE(1) = y0;
MASS(1) = m_dry + m_prop + m_motor;
INT_THRUST(1) = T0;

t = 0;
landInd = 0;
burnoutInd = 0;
ejectInd = 0;
for k = 2:numSteps
    t = t + dt;
    TIME(k) = t;
    %Get the thrust at that time
    if(t <= THRUST(size(THRUST, 1)))
        T = getThrust(t, THRUST);
        MASS(k) = MASS(k - 1) - m_dot_avg*dt;
        burnoutInd = k;
    else
        T = 0;
        MASS(k) = MASS(k - 1);
    end
    
    %Get the physical constants
    rho = getDensity(ALTITUDE(k - 1)); % [kg/m^3]
    g = getGravity(ALTITUDE(k - 1), g0); % [m/s^2]
    
    %Force Analysis
    
    %Drag
    if(t >= THRUST(size(THRUST, 1)) + delay) % Parachute opens
        D = .5*CD_para*rho*VELOCITY(k - 1).^2*A_parachute; % [N]
    else %Parachute has not opened yet
        D = .5*CD*rho*VELOCITY(k - 1).^2*A_body; % [N]
        ejectInd = k;
    end
    
    %Check Direction of Travel
    if(VELOCITY(k - 1) <= 0)
        F_net = T + D - MASS(k)*g*cos(theta); % [N]
    else
        F_net = T - D - MASS(k)*g*cos(theta); % [N]
    end
    
    %Get the acceleration
    a_z = F_net/MASS(k); % [m/s^2]
    
    %Add to vectors
    F_NET(k) = F_net;
    if(F_net <= 0 && t < .5)
        ALTITUDE(k) = 0;
        VELOCITY(k) = 0;
        ACCELERATION(k) = 0;
    elseif(landInd ~= 0)
        ALTITUDE(k) = 0;
        VELOCITY(k) = 0;
        ACCELERATION(k) = 0;
    else
        ACCELERATION(k) = a_z;
        VELOCITY(k) = VELOCITY(k - 1) + a_z*dt;
        ALTITUDE(k) = ALTITUDE(k - 1) + VELOCITY(k)*dt;
        if(ALTITUDE(k) < 0)
            landInd = k;
        end
    end
end

%Analysis
[apo, apo_index] = max(ALTITUDE);
[v_max, vmax_index] = max(VELOCITY);
[a_max, amax_index] = max(ACCELERATION);
y = linspace(-500,800);

off_rail_ind = 1;
while(ALTITUDE(off_rail_ind)*3.28 <= rail_height)
    off_rail_ind = off_rail_ind + 1;
end

%Plot the Results
figure(1) %Altitude
plot(TIME, ALTITUDE*3.28);
hold on;
plot(THRUST(size(THRUST, 1)), ALTITUDE(round(THRUST(size(THRUST, 1),1)/dt + 1))*3.28, 'b.'); % Burnout
plot(TIME(apo_index), apo*3.28, 'g*'); % Apogee
plot(delay + THRUST(size(THRUST, 1)),... % Ejection
    ALTITUDE(round((delay + THRUST(size(THRUST, 1),1))/dt))*3.28, 'r*');
title(['Altitude Plot: C_D = ' ...
    num2str(CD) ', C_D Parachute = ' num2str(CD_para)]);
xlabel('Time [s]');
ylabel('Altitude [ft]');
legend('Simulated Altitude', 'Burnout', 'Apogee', 'Ejection Charge');
xlim([0 TIME(landInd)]); ylim([0 1.1*apo*3.28]);
text(TIME(apo_index) + TIME(landInd)/20, apo*3.28, ['Apogee: ' num2str(apo*3.28) ' ft']);
grid on;

figure(2) %Velocity
plot(TIME, VELOCITY*3.28);
hold on;
plot(THRUST(size(THRUST,1),1)*ones(1, 100), y, 'r-'); %Burnout
plot(TIME(apo_index)*ones(1, 100), y, 'b-'); % Apogee
plot((THRUST(size(THRUST,1),1) + delay)*ones(1, 100), y, 'g-'); % Ejection
title(['Velocity Plot: C_D = ' ...
    num2str(CD) ', C_D Parachute = ' num2str(CD_para)]);
xlabel('Time [s]');
ylabel('Velocity [ft/s]');
if(CD_para ~= 0)
    xlim([0 ((THRUST(size(THRUST,1),1) + delay) + 2)])
    text(5, VELOCITY(landInd),[...
    'Touchdown Velocity: ' num2str(VELOCITY(landInd)) ' ft/s']);
else
    xlim([0 TIME(landInd - 3)]);
    text(TIME(ceil(landInd/2)), VELOCITY(landInd),[...
    'Touchdown Velocity: ' num2str(VELOCITY(landInd)) ' ft/s']);
end
ylim([1.1*min(VELOCITY)*3.28 1.1*max(VELOCITY)*3.28]);
text(TIME(vmax_index), v_max*3.28, ['Max Velocity: ' num2str(v_max*3.28) ' ft/s']);
text(TIME(off_rail_ind), VELOCITY(off_rail_ind)*3.28,...
    ['Off the Rail Velocity: ' num2str(VELOCITY(off_rail_ind)*3.28) ' ft/s']);
legend('Velocity', 'Burnout', 'Apogee', 'Ejection Charge');
grid on;

figure(3) %Acceleration
plot(TIME, ACCELERATION*3.28);
hold on;
plot(THRUST(size(THRUST,1),1)*ones(1, 100), y, 'r-'); %Burnout
plot(TIME(apo_index)*ones(1, 100), y, 'b-'); % Apogee
plot((THRUST(size(THRUST,1),1) + delay)*ones(1, 100), y, 'g-'); % Ejection
title(['Acceleration Plot: C_D = ' ...
    num2str(CD) ', C_D Parachute = ' num2str(CD_para)]);
xlabel('Time [s]');
ylabel('Acceleration[ft/s^2]');
xlim([0 ((THRUST(size(THRUST,1),1) + delay) + 2)]);
ylim([1.1*min(ACCELERATION)*3.28 1.1*max(ACCELERATION)*3.28]);
legend('Acceleration', 'Burnout', 'Apogee', 'Ejection Charge');
text(TIME(amax_index), a_max*3.28, ['Max Acceleration: ' num2str(a_max*3.28) ' ft/s^2']);
grid on;

%Print useful data to screen
fprintf('Apogee: %d ft\n', apo*3.28);
fprintf('Flight Duration: %d s\n', TIME(landInd));
fprintf('Burnout Velocity: %d ft/s\n', VELOCITY(burnoutInd + 1)*3.28);
fprintf('Coast Change in Alt: %d ft\n',...
    (ALTITUDE(ejectInd + 1) - ALTITUDE(burnoutInd - 1))*3.28);