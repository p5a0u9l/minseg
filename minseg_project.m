%% MinSeg Project - EE 547 (PMP) - Winter 2015
% prepared by Paul Adams
%

%% Initialization
close all
digits(3);
set(0, 'defaultTextInterpreter', 'latex'); 
format shortG
numerical_precision = 1e-9;
syms a x I_p m_p L r_w I_w m_w r_w
syms s k_t R V k_b 

%% 4.1 Dynamical Model of the MinSeg Robot
%%
% <html> <h3> Physical Parameters. </h3> </html>
g = 9.80665;    % [m/s^2]
k_t = 0.3233;   % [Nm/a]
k_b = 0.4953;   % [Vs/rad]
R = 5.2628;     % [Ohms]
L = 0.05;       % [m]   - guess
m_p = 0.3;      % [kg]  - guess
m_w = 0.1;      % [kg]  - guess
r_w = 0.01;     % [m]   - guess
% assume a filled circular area
I_w = pi/4*r_w^4; % [m^4]   - http://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
w_arduino = 0.05363; % [m] - width of arduino, http://www.adafruit.com/product/191
h_arduino = 0.01529; % [m] - height of arduino, http://www.adafruit.com/product/191
l_arduino = 0.10198; % [m] - length of arduino, http://www.adafruit.com/product/191
% assuming a filled rectangular area
I_p = w_arduino*l_arduino^3/12; % [m^4] - http://en.wikipedia.org/wiki/List_of_area_moments_of_inertia

%%
% <html> <h3> State-space Matrices. </h3> </html>
Arow12 = (g*L*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow22 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow24 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*r_w*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow41 = (g*L^2*m_p^2*r_w^2)/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow42 = -k_b*k_t*r_w*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow44 = -k_b*k_t*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A = [0, 1, 0, 0; Arow12, Arow22, 0, Arow24; 0, 0, 0, 1; Arow41, Arow42, 0, Arow44];    
render_latex(['A = ' latex(vpa(A, 3))], 10, 1)
%% 
Brow2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow3 = -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B = [0; Brow2; 0; Brow3];
render_latex(['B = ' latex(vpa(B, 3))], 10, 1)
n = size(A, 1);
C = eye(n);
D = zeros(n, 1);
%%
% <html> <h3> Transfer Function. </h3> </html>
sys = ss(A, B, C, D);
[num, den] = ss2tf(A, B, C, D);
for i = 1:n
    G(i, :) = vpa(poly2sym(num(i, :), s), 2)/vpa(poly2sym(den, s), 2);
end
render_latex(['\hat{G}(s) = ' latex(vpa(G, 2))], 12, 1.5)
%%
% <html> <h3> Characteristic Polynomial and eigenvalues. </h3> </html>
Delta = vpa(charpoly(A, s), 2);
render_latex(['\Delta(\lambda) = ' latex(vpa(Delta, 2))], 12, 0.5)
lambda = eig(A)
%%
% <html> <h3> Check if the system is asymptotically stable. </h3> </html>
if all(real(eig(A)) < 0)
    disp('System is Asymptotically stable')
else
    disp('System is Not Asymptotically stable')
end
%%
% <html> <h3> Step Response of open-loop system. </h3> </html>
[y, t, x] = step(sys);
f = figure;
f.Position(3) = 1.5*f.Position(3);
plot(t, y)
xlabel('time [s]')
title('Step-input response of open-loop system')
legend('\alpha', '\dot{\alpha}', 'x', '\dot{x}', 'Location', 'northwest')
%%
% <html> <h3> Check if the system is controllable by the rank of controllability matrix by MATLAB <i>rank</i> function.</h3> </html>
Cm = ctrb(sys.a, sys.b);
if rank(Cm) >= n
    disp('System is controllable')
else
    disp('System is not controllable')
end
ccf = canon(sys);

%%
close all
