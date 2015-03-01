%% MinSeg Project - EE 547 (PMP) - Winter 2015
% prepared by Paul Adams
%

%% Initialization
addpath('simulink')
close all
digits(3);
set(0, 'defaultTextInterpreter', 'latex'); 
format shortG
numerical_precision = 1e-9;
syms a x I_p m_p L r_w I_w m_w r_w
syms s k_t R V k_b 

%% 4.1 Dynamical Model of the MinSeg Robot
%%
% <html> <h3> Step 2 Physical Parameters. </h3> </html>
g = 9.80665;    % [m/s^2]
k_t = 0.3233;   % [Nm/a]
k_b = 0.4953;   % [Vs/rad]
R = 5.2628;     % [Ohms]
L = 0.11;       % [m]   - demonstrated balance point with 6 AA batteries
m_p = 0.3;      % [kg]  - guess
m_w = 0.1;      % [kg]  - guess
r_w = 0.016;    % [m]   - 5/8", measured
% assume a filled circular area (x2 for inertia of both wheels)
I_w = m_w*r_w^2/2;   % [kg-m^2]   - http://en.wikipedia.org/wiki/List_of_moments_of_inertia
w_arduino = 0.05363; % [m] - width of arduino, http://www.adafruit.com/product/191
h_arduino = 0.01529; % [m] - height of arduino, http://www.adafruit.com/product/191
h_p = 0.2;           % [m] - height of pendulum, from top of PCB to wheel axis, measured
l_arduino = 0.10198; % [m] - length of arduino, http://www.adafruit.com/product/191
% assuming a filled rectangular area
%I_p = w_arduino*l_arduino^3/12; % [m^4] - http://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
% assuming rod length L and mass m
%I_p = m_p * h_p^2 / 3; % [kg-m^2] - http://en.wikipedia.org/wiki/List_of_moments_of_inertia
% assuming point mass
I_p = m_p * L^2; % [kg-m^2] - http://en.wikipedia.org/wiki/Moment_of_inertia

%%
% <html> <h3> Step 1 State-space Matrices. </h3> </html>
Arow12 = (g*L*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow22 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow24 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*r_w*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow41 = (g*L^2*m_p^2*r_w^2)/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow42 = -k_b*k_t*r_w*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow44 = -k_b*k_t*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A = [0, 1, 0, 0; Arow12, Arow22, 0, Arow24; 0, 0, 0, 1; Arow41, Arow42, 0, Arow44];    
%render_latex(['A = ' latex(vpa(A, 3))], 10, 1)

%% 
Brow2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow3 = -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B = [0; Brow2; 0; Brow3];
%render_latex(['B = ' latex(vpa(B, 3))], 10, 1)
n = size(A, 1);
C = eye(n);
D = zeros(n, 1);

%%
% <html> <h3> Step 3 Transfer Function. </h3> </html>
sys = ss(A, B, C, D);
[num, den] = ss2tf(A, B, C, D);
for i = 1:n
    G(i, :) = vpa(poly2sym(num(i, :), s), 2)/vpa(poly2sym(den, s), 2);
end
%render_latex(['\hat{G}(s) = ' latex(vpa(G, 2))], 12, 1.5)
%%
% <html> <h3> Step 4 Characteristic Polynomial and eigenvalues. </h3> </html>
Delta = vpa(charpoly(A, s), 2);
%render_latex(['\Delta(\lambda) = ' latex(vpa(Delta, 2))], 12, 0.5)
lambda = eig(A)

%%
% <html> <h3> Step 5 Check if the system is asymptotically stable. </h3> </html>
if all(real(eig(A)) < 0)
    disp('System is Asymptotically stable')
else
    disp('System is Not Asymptotically stable')
end
%%
% <html> <h3> Step 6 Find the poles of the transfer function </h3> </html>
[~,poles_minseg,~]=zpkdata(sys) %todo - BIBO stable?

%% 4.2 Controllability and Observability of the System
% <html> <h3> Step 7 Check if the system is controllable by the rank of controllability matrix by MATLAB <i>rank</i> function.</h3> </html>
Cm = ctrb(sys.a, sys.b);
if rank(Cm) >= n
    disp('System is controllable')
else
    disp('System is not controllable')
end

%%
% <html> <h3> Step 8 Analyze the observability of the linearized system.</h3> </html>
Om = obsv(sys.a, sys.c);
if rank(Cm) >= n
    disp('System is observable')
else
    disp('System is not observable')
end

%%
% <html> <h3> Step 9 Transform the linearized system into a _controllable canonical form_ and _observable canonical form_.</h3> </html>
ccf = canon(sys,'companion'); %todo observable canonical form

%% 4.3 State Estimator
% <html> <h3> Step 10 Develop a closed loop state estimator for the open loop system.</h3> </html>
poles_obsv = -6*abs(poles_minseg{3}); %todo fix this
L = place(transpose(A),transpose(C),poles_obsv);

%%
% <html> <h3> Step 11 Develop a Simulink model of the linearized system.</h3> </html>
xini = [0 0 0 0];
xhatini = [0 0 0 0];
sim('step_11');
figure
subplot(3,1,1)
plot(time,x)
subplot(3,1,2)
plot(time,xhat)
subplot(3,1,3)
plot(time,y)

%%
% <html> <h3> Step Response of open-loop system. </h3> </html>
[y, t, x] = step(sys);
f = figure;
%f.Position(3) = 1.5*f.Position(3);
plot(t, y)
xlabel('time [s]')
title('Step-input response of open-loop system')
legend('\alpha', '\dot{\alpha}', 'x', '\dot{x}', 'Location', 'northwest')


%%
%close all
