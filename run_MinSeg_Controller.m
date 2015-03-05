clear
close
addpath('simulink')

%% Physical Parameters
g = 9.80665;    % [m/s^2]
k_t = 0.3233;   % [Nm/a]
k_b = 0.4953;   % [Vs/rad]
R = 5.2628;     % [Ohms]
L = 0.11;       % [m]   - demonstrated balance point with 6 AA batteries
m_brick_bat = 0.249; % measured in class
m_wheel = 0.014; % measured in class
m_motor = 117 - 2*m_wheel; 
m_p = m_brick_bat + m_motor;      % [kg]  - guess
m_w = 2*m_wheel;      % [kg]  - guess
% m_p = 0.3;      % [kg]  - guess
% m_w = 0.1;      % [kg]  - guess
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

%% Open-Loop State-space Matrices
Arow12 = (g*L*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow22 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow24 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*r_w*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow41 = (g*L^2*m_p^2*r_w^2)/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow42 = -k_b*k_t*r_w*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow44 = -k_b*k_t*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow3 = -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));

A = [0, 1, 0, 0; Arow12, Arow22, 0, Arow24; 0, 0, 0, 1; Arow41, Arow42, 0, Arow44];    
B = [0; Brow2; 0; Brow3];
N = size(A, 1);
C = eye(N);
D = zeros(N, 1);
sys = ss(A, B, C, D);
[num, den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Cm = ctrb(sys.a, sys.b);

%% Control-Canonical Form 
alpha = den(2:end); % denominator coefficients of G(s)
Cm_bar_inv = [1, alpha(1), alpha(2), alpha(3); 
              0, 1, alpha(1), alpha(2); 
              0, 0, 1, alpha(1);
              0, 0, 0, 1];
Q = Cm*Cm_bar_inv; 
A = round(Q\sys.A*Q*1e5)/1e5; % remove approximation errors close to zero
B = round(Q\sys.B*1e5)/1e5;
C = round(sys.C*Q*1e5)/1e5;
D = D;
ccf = ss(A, B, C, D);

%% Select a state-space
selection = 'ccf';
switch selection
    case 'ccf'
        A = ccf.a;
        B = ccf.b;
        C = ccf.c;
        D = ccf.d;
    case 'sys'
        A = sys.a;
        B = sys.b;
        C = sys.c;
        D = sys.d;
end

%% Proportional Feedback Controller
poles_fbkCtrl = [-10+5j, -10-5j, -12+1j,-12-1j];  %PRA try picking arbitray e-vals?
K = place(A, B, poles_fbkCtrl);
A = A - B*K;

%% simulation 
cd simulink
xini = [0 0 0 0];
sim_sample_freq = 0.3;
% sim('MinSeg_Controller', inf)
set_param('MinSeg_Controller', 'SimulationCommand', 'start')