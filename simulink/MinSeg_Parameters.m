%% constants or conversion factors:
R = 5.262773292;     % ohms
k_b = 0.4952900056;    % Vs/rad
k_t = 0.3233728703;    % Nm/A
Bm = 0.0006001689451; % Nms/rad (viscous friction)
Lm = 0.0047;          % H    
Jm = 0.001321184025;  % kgm^2 (combined J)
Tf = 0.007299397206;  % Nm (coulomb friction)
encoder_counts=720;   % number of counts (if using quad encoding)
RPM_MAX = 170;        % spec sheet max RPM
RADSEC2RPM = 60/(2*pi);       % radians/sec to RPM
RAD2C=encoder_counts/(2*pi);  % conversion from radians to counts
C2RAD=1/RAD2C; 
C2DEG=360/encoder_counts;
R2D=180/pi;
D2R=1/R2D;
g = 9.80665;    % [m/s^2]
GyroS = 1/112;  % gyro scaling factor

%% MinSeg Parameters
r_w = 0.016;    % [m] - 5/8", measured
m_board_bat = 0.250; % measured in class
m_board = 0.120;
m_wheel = 0.017; % measured in class
m_motor = 0.084; 
m_battery = 6/0.252; 
num_bats = 3;
m_p = m_board_bat + m_motor + num_bats*m_battery;      % [kg]
switch num_bats
    case 0
        L_com = 0.095;       % [m]   - demonstrated balance point with 6 AA batteries
        Vsupply = 5;  % 9 volts
    case 3
        L_com = 0.1;
        Vsupply = 5;  % 9 volts
    case 6
        L_com = 0.11;       % [m]   - demonstrated balance point with 6 AA batteries
        Vsupply = 9;  % 9 volts
end

%% Dependant Parameters
DCB2V = Vsupply/255;   % Duty cycle bits to PWM (volts)
V2DCB = 1/DCB2V;       % Volts to duty cycle in bits.m_w = 2*m_wheel;      % [kg]
I_w = m_w*r_w^2/2;   % [kg-m^2]   - http://en.wikipedia.org/wiki/List_of_moments_of_inertia
% assuming point mass
I_p = m_p * L_com^2; % [kg-m^2] - http://en.wikipedia.org/wiki/Moment_of_inertia

%% Open-Loop State-space Matrices
Arow12 = (g*L_com*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow22 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L_com + r_w)))/(R*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow24 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L_com + r_w)))/(R*r_w*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow41 = (g*L_com^2*m_p^2*r_w^2)/(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
Arow42 = -k_b*k_t*r_w*(I_p + L_com*m_p*(L_com + r_w))/(R*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Arow44 = -k_b*k_t*(I_p + L_com*m_p*(L_com + r_w))/(R*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L_com + r_w))))/(R*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
Brow3 = -(k_t*r_w*(I_p+ L_com*m_p*(L_com + r_w)))/(R*(I_w*(I_p + L_com^2*m_p) + (L_com^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));

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
selection = 'sys';
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

%% LQR Feedback Controller
Q = C'*C;
% Q(1, 1) = 200;
% Q(3, 3) = 200;
R = 1; 
KLQR = lqr(A, B, Q, R);
% KLQRC=KLQR;
% KLQRC = KLQR.*[-r_w -r_w 1 1]; % combined LQR Gain with wheel radius scaling
% Ki=1*[-r_w -r_w*0 1 1*0];

TS = .03;

% combine some constants:
GS = GyroS*D2R;
ES = -C2DEG*D2R;

% automatic calibration of gyro offset
% discrete 1st order filter recurrence relation - 
% discrete-time implementation of a low-pass filter is 
% the exponentially-weighted moving average
% alpha=.5 (time contant is equal to sampling period)
% alpha<.5 (time constant is larger than sampling period) tau ~ TS/alpha

% .006 requires 800 samples to get to steardy state (see gyro filter design)
% 1st order system take 3Tau to get to .95 ss value, 3.9Tau, 98, 4.56tau 99
% tau ~ .0025/.006 = .4167 = 3*.4167 = 1.25 to get to .95  - we need .99
% tau ~ .0025/.006 = .4167 = .4167*.456 = 1.9sec to get to 99% ss
a_go = .006;  % alpha for initial gyro offset calibration (2 sec for ss)

% Start time for balancing - Gyro calibartion time:
tstart = 2;
K = KLQR;

%% MPU5060
return
% Low pass filter for output:
a_ov = .7;  % on USB
a_ov= .5; % on battery
a_ov=1;  % no filter.

