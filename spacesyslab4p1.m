%%Setup
% A script written in MATLAB with the objective of describing certain
% characteristics of the entry of the Mars Pathfinder probe. 
% TODO: Convert to python (in progress)
% TODO: Check with professor about part H - I would be shocked if I 
%       did that part properly.
% Last Edit: 29 October 2016

% Useful Constants and Quantities:
gE  = 9.81 ;                     % Acc. due to gravity, Earth, m/s2
gM  = 3.69 ;                     % Acc. due to gravity, Mars, m/s2
RM  = 3397 * 1000 ;              % Radius of Mars, m
sig = 5.67e-8 ;                  % Stefan-Boltzmann Constant, W/m2*K4
k   = 1.9027e-4 ;                % Sutton-Graves const. for Mars
R_univ = 8314.472 ;              % Universal Gas Constant, J/(K*kmol)

% Problem 1 Givens
d = 2.65 ;                       % Diameter of entry capsule, m
delta = 70 ;                     % Angle of sphere-cone heatshield, degrees
CD = 1.65 ;                      % Hypersonic Drag Coefficient
m = 603 ;                        % Entry mass, kg

% Problem 1 Assumptions
rn = 0.665 ;                     % Effective nose radius, m
rho0 = 0.0525 ;                  % Adjusted density at 'sea' level, kg/m3
T = 150 ;                        % Isothermal Temp, K (const. throughout)
M = 43.27 ;                      % Mean molecular weight, kg/kmol

% Entry Initial conditions
V_atm = 4.9 * 1000 ;             % m/s
h_atm = 125 * 1000 ;             % m
gam_atm = 17.0 ;                 % degrees
                                 % Heatshield material: SLA-561
T0_TPS = 25 + 273.15 ;           % Init. uniform temp distribution in the 
                                 % Thermal Protection System, K                         
T_max = 200 + 273.15 ;           % Maximum allowable temp of bondline, K

%% Part A
% Find the ballistic coefficient and the Allen-Eggers constant. 

% beta = m/(CD * A)
% By EDL and TPS Slide 9, for reference area we use the projected area of
% the cone. 

A = pi * (d / 2) ^ 2 ;           % CS Area of capsule
beta = m / ( CD * A ) ;          % Ballistic coefficient
fprintf('\nThe ballistic coefficient is : %f \n',beta);

% Allen-Eggers Constant C = - rho0 * H / (2 * beta * sin(gam_atm))
% Need scale height H = R * T0 / (g0)
R = R_univ / M ;                               % Spec. Gas Const., J / kgK
H = R * T / gM ;                               % Atmo. scale height, m
C = - rho0 * H / (2 * beta * sind(gam_atm)) ;  % A-E Const
fprintf('Atmospheric Scale Height: %f m \n',H) ;
fprintf('Allen-Eggers Constant: %f \n', C) ;

%% Part B - Allen-Eggers Velocity Profile
% Use Allen-Eggers eq. to develop a velocity profile as function of
% altitude. Plot results with v on x, h on y, km/s and km respectively. 
h = 0:1:h_atm ;    % Array for h in m w/ step size in m. 
V = V_atm .* exp( C .* exp(-h ./ H)) ;
plot(V ./1000 , h ./ 1000)
hold on
set(gca,'XDir','Reverse') ;          % I figured it'd be neat to plot the 
                                     % plot the descent from left to right. 
xlabel('Velocity (km/s)') ;
ylabel('Altitude (km)');
title('Altitude vs Velocity of Pathfinder EDL via Allen-Eggers Model') ;

%% Part C - Maximum Deceleration
% Find the altitude of max. deceleration. Express the magnitude of max.
% deceleration in Earth gs. 

h_nMAX = H * log(-2 * C) ;                      % Alt. of max. decel., m
nMAX = V_atm .^ 2 * sind(gam_atm) / (53.3 * H) ;% Mag peak decel, Earth gs.
fprintf('Altitude of max deceleration: %f km \n', h_nMAX / 1000) ;
fprintf('Max deceleration: %f Earth gs \n', nMAX)

%% Part D - Equations of Motion
% Integrate EOM using ode45. Plot velocity profile against profile found in
% (b). Comment on the quality of the assumptions made going from EOM to
% Allen-Eggers solution. Estimate time from the Phoenix EDL profile shown
% in lecture EDL and TPS, slide 4.
% dV/dt = -(rho * V^2)/(2 * beta) + g * sin(gamma)
% dh/dt = -V * sin(gamma)
% dgamma/dt = - (V * cos(gamma)) / (Rp + h) + g * cos(gamma) / V

tf = 221 ;                               % Time of chute deployment, s.
TOF = [0:1:tf] ;                         % Time of flight, s.

VHG = [ V_atm h_atm gam_atm ]' ;         % Initial Conditions

% The anonymous function that I'll be calling in the ode45() function. I
% realize it's not the most readable, but it works. Each row corresponds to
% the differential equation for dV/dt, dh/dt and dgamma/dt respectively.
% For dV/dt, I had to define rho according to the exponential atmosphere
% model within the anonymous function in order for it to not be assumed
% constant. Also, I chose to remove the whitespace from the first and
% last rows to ensure that they fit within the confines of what the MATLAB
% editor considers the proper length for a line. I realize that I choose 
% the strangest things to be particular about. 
f=@(t,VHG)[-((rho0*exp(-VHG(2)/H))*VHG(1).^2)/(2*beta)+gM.*sind(VHG(3));...
           -VHG(1) .* sind(VHG(3))                                     ;...
           -(VHG(1).*cosd(VHG(3)))./(RM+VHG(2))+gM.*cosd(VHG(3))./VHG(1)];

[t,x] = ode45(f,TOF,VHG);

figure

% Get the velocity according to the Allen-Eggers equation using the heights
% found from solving the differential equations. 
V_ae = V_atm .* exp( C .* exp(-x(:, 2) / H));
plot(x(:,1) ./1000, x(:,2) ./ 1000)
hold on
plot(V_ae ./ 1000, x(:,2) ./ 1000)
ylabel('Altitude (km)');
xlabel('Velocity (km/s)');
title('Comparison of EOM Solutions to Allen-Eggers Model');
legend('EOM','Allen-Eggers','Location','NorthEast');
set(gca,'XDir','Reverse') ;         % Again set x axes to reverse so that 
                                    % the plot is indicative of what's
                                    % happening with increasing time.
                                    
%% Part E - Heat Rate vs Altitude
% For the integrated and Allen-Eggers solution, use the Sutton-Graves
% equation to estimate the stagnotion point heat rate as a function of
% altitude. Plot the altitude (km) vs the heat rate (W/cm2)
% The heat rate is given in W/m2 by qdot = sqrt(k(rho/rn)) * V^3

% By slide 29, approximate effective nose radius as rn/2.
% Heat rate in W/m2 for Allen-Eggers solutions
% Need to recompute A-E solutions for heights found by integration.
V_re = V_atm .* exp( C .* exp(-x(:,2) ./ H)) ;
qd_AE = k .* sqrt(((rho0 .* exp(-x(:, 2) ./ H)) ./ (rn/2))) .* V_re .^ 3;

% Heat rate in W/m2 for integrated equations
qd_int=k .* sqrt(((rho0 .* exp(-x(:, 2) ./ H)) ./ (rn/2))) .* x(:,1) .^ 3 ;

figure
plot(qd_int .* 1e-4, x(:, 2) ./ 1000, qd_AE .* 1e-4, x(:, 2) ./1000)
xlabel('Heating Rate (W/cm^2)'),  ylabel('Altitude (km)')
title('Heat rate via Sutton-Graves vs Altitude')
legend('Integrated Solution','Allen-Eggers')

%% Part F - Peak Heat Rate and Heat Load
% Use analytic expressions to estimate the peak stagnation point heat rate
% and the integrated heat load. Compare the peak heat rate estimate to the
% maximum values from part E.

% Peak heat rate, assuming constant flightpath angle. 
qd_max=k * sqrt((beta * sind(gam_atm) / (3 * (rn/2) * H))) * ...
       (V_atm ^3 / exp(.5)) ;

% Total Integrated Heat Load Estimate, assuming constant flightpath angle. 
hl_est = k * V_atm ^ 2 * sqrt(pi * H * beta / ((rn/2) * sind(gam_atm))) ;

% Get peak heat rates from part E and compare them to the above peak heat
% rate estimate. Compare both in terms of simple magnitude and percent
% difference, treating the values from part E as the 'true' values. 
qd_AE_max = max(qd_AE) ; 
qd_int_max = max(qd_int) ;

diffAE = abs(qd_max - qd_AE_max) ;
diffINT = abs(qd_max - qd_int_max) ;
pDiffAE = diffAE / ( 0.5 * ( qd_max + qd_AE_max)) * 100 ;
pDiffINT = diffINT / ( 0.5 * ( qd_max + qd_int_max)) * 100 ;
fprintf('Peak Heat Rate with constant flightpath angle: %f\n',qd_max) ;
fprintf('Peak Heat Rate of integrated solution: %f\n',qd_int_max) ;
fprintf('Peak Heat Rate of Allen-Eggers solution: %f\n',qd_AE_max) ;
fprintf('INTEGRATED: Diff = %f, Percent Diff = %f\n',diffINT,pDiffINT) ;
fprintf('ALLEN-EGGERS: Diff = %f, Percent Diff = %f\n',diffAE,pDiffAE) ;

%% Part G - Temperature
% Assume outer surface is in radiative equilibrium. What is the TPS outer
% surface temp as a function of time? Use the integrate results to get
% this. Plot with time on the x-axis, temp on y-axis in Celsius. 
Eps = 0.9 ;                     % Emissivity of SLA-561
Temp = (qd_int ./ (sig .* Eps)) .^ (1 / 4) ;
Temp = Temp + 273.15 ;                % Convert to Celsius

figure
plot(t,Temp), xlabel('Time (s)'), ylabel('Temperature (deg C)')
title('Temperature vs Time')

%% Part H - Mass of TPS
% Estimate the mass of the TPS material vs the integrated heat load. What
% is the TPS material mass? What is the mass of the underlying TPS
% structure? What is the total TPS mass? What is the TPS thickness (assume
% a 70-degree cone)?
% Use the heat load estimate from part F and the equation from the plot on
% slide 34 of EDL and TPS lecture.

MF = (0.091 * (hl_est * 1e-4) ^ 0.51575) / 100 ;  % Mass fraction in %
                                                  % J/m2 ---> J/cm2
MF_struc = [ 0.1 0.15] ;               % Mass fraction of support structure
fprintf('The mass of the TPS material is %f kg\n', MF * m)
fprintf('The mass of the TPS structure is between %f kg and %f kg\n', ...
        MF_struc(1) * m, MF_struc(2) *m)

rhoTPS = 0.288 * 1000 ;             % Density of TPS material, kg/m3
ATPS = pi * d ^ 2 / 4 ;             % Area of heatshield
thickness = (MF * m) / (rhoTPS * ATPS) ;
fprintf('The thickness of the TPS is: %f cm\n', thickness * 100)

%% Part I - Acceleration
% Check the baseline peak deceleration of 9.2 Earth gs against Allen-Eggers
% and integrated results.
n_base = 9.2 ;
dVdtINT = - ((rho0 .* exp(-x(:, 2) / H)) .* x(:,1) .^2) ./ (2 .* beta) ...
          + gM * sind(x(:, 3)) ;
n_max_int = max(abs(dVdtINT)) / gE ;
fprintf('The difference between the Allen-Eggers and baseline peak deceleration is: %f gs\n', nMAX - n_base);
fprintf('The difference between the integrated and baseline peak deceleration is: %f gs\n', n_max_int - n_base);

%% Part J - Flightpath Angle
% Plot the flightpath angle in degrees vs the altitude in km. Discuss the
% cause of errors in A-E heating estimates and the max deceleration
% estimates.

figure
plot(x(:,2) ./ 1000, x(:,3))
xlabel('Altitude (km)') ;
ylabel('Flightpath Angle (deg)');
title('Flightpath Angle vs Altitude');
set(gca,'XDir','Reverse') ;   

%% Bonus - Distance traveled 
tf = 221 ;                               % Time of chute deployment, s.
TOF = [0:1:tf] ;                         % Time of flight, s.
% Create new initial conditions for x0
VHG = [ V_atm h_atm gam_atm 0.0]' ;         % Initial Conditions


f=@(t,VHG)[-((rho0*exp(-VHG(2)/H))*VHG(1).^2)/(2*beta)+gM.*sind(VHG(3));...
           -VHG(1) .* sind(VHG(3))                                     ;...
           -(VHG(1).*cosd(VHG(3)))./(RM+VHG(2))+gM.*cosd(VHG(3))./VHG(1);...
           VHG(1) * cosd(VHG(3))                                        ];
[t,x] = ode45(f,TOF,VHG);
figure
plot(x(:,4) ./ 1000, x(:,2) ./1000)
xlabel('Distance Traveled Downrange (km)')
ylabel('Altitude (km)')
title('Altitude vs Distance Traveled')
